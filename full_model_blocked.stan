data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> n_genus;
    int<lower=0> n_B;
    int<lower=0> n_B_g;
    int<lower=0> n_B_T;
    int<lower=0> n_period;
    int<lower=0> n_blocks;
    real<lower=0> sigma_obs_lower;
    int<lower=0> tree_ID[n_tree];
    int<lower=0> plot_ID[n_tree];
    int<lower=0> site_ID[n_tree];
    int<lower=0> genus_ID[n_tree];
    int<lower=0> t0[n_tree];
    int<lower=0> tf[n_tree];
    int<lower=0> n_miss;
    int<lower=0> n_obs;
    int<lower=0> obs_indices_tree[n_obs];
    int<lower=0> obs_indices_period[n_obs];
    int<lower=0> miss_indices_tree[n_miss];
    int<lower=0> miss_indices_period[n_miss];
    int<lower=0> block_n_periods[n_blocks];
    int<lower=0> block_start_row[n_blocks];
    int<lower=0> block_end_row[n_blocks];
    vector[n_obs] dbh_obs;
    vector[n_tree] WD;
    vector[n_tree] WD_sq;
    vector[n_tree] elev;
    matrix[n_tree, n_period + 1] temp;
    matrix[n_tree, n_period + 1] precip;
    matrix[n_tree, n_period + 1] precip_sq;
}

parameters {
    matrix[n_tree, n_period + 1] dbh_latent;
    vector[n_miss] dbh_miss;
    vector[n_B] B;
    matrix[n_site, n_B_T] B_T;
    matrix[n_B_g, n_genus] B_g_std;
    vector[n_B_g] gamma_B_g;
    cholesky_factor_corr[n_B_g] L_rho_B_g;
    vector<lower=0>[n_B_g] sigma_B_g_sigma;
    vector[n_tree] int_ijk_std;
    vector[n_plot] int_jk_std;
    vector[n_site] int_k_std;
    real<lower=sigma_obs_lower> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_int_ijk;
    real<lower=0> sigma_int_jk;
    real<lower=0> sigma_int_k;
}

transformed parameters {
    matrix[n_tree, n_period + 1] dbh;
    matrix[n_tree, n_period + 1] dbh_latent_sq;
    vector[n_tree] int_ijk;
    vector[n_plot] int_jk;
    vector[n_site] int_k;
    matrix[n_genus, n_B_g] B_g;

    // Handle missing data
    for (n in 1:n_miss) {
        dbh[miss_indices_tree[n], miss_indices_period[n]] <- dbh_miss[n];
    }

    for (n in 1:n_obs) {
        dbh[obs_indices_tree[n], obs_indices_period[n]] <- dbh_obs[n];
    }

    // Matt trick- see http://bit.ly/1qz4NC6
    int_ijk <- sigma_int_ijk * int_ijk_std; // int_ijk ~ normal(0, sigma_int_ijk)
    int_jk <- sigma_int_jk * int_jk_std; // int_jk ~ normal(0, sigma_int_jk)
    int_k <- sigma_int_k * int_k_std; // int_k ~ normal(0, sigma_int_k)

    B_g <- transpose(rep_matrix(gamma_B_g, n_genus) + diag_pre_multiply(sigma_B_g_sigma, L_rho_B_g) * B_g_std);

    // TODO: This could be made faster by not squaring cells that are no data
    for (i in 1:n_tree) {
        for (j in 1:(n_period+1)) {
            dbh_latent_sq[i][j] <- square(dbh_latent[i][j]);
        }
    }
}

model {
    int n_rows;

    //########################################################################
    // Fixed effects
    B ~ normal(0, 10);

    //########################################################################
    // Temperature model (per site)
    to_vector(B_T) ~ normal(0, 10);

    //########################################################################
    // Observation and process error
    sigma_obs ~ cauchy(0, 1);
    sigma_proc ~ cauchy(0, 1);

    //########################################################################
    // Nested random effects
    int_ijk_std ~ normal(0, 1); // Matt trick
    sigma_int_ijk ~ cauchy(0, 1);

    int_jk_std ~ normal(0, 1); // Matt trick
    sigma_int_jk ~ cauchy(0, 1);

    sigma_int_k ~ cauchy(0, 1);
    int_k_std ~ normal(0, 1); // Matt trick

    //########################################################################
    // Correlated random effects at genus level (crossed). Based on Gelman and 
    // Hill pg. 377-378, and http://bit.ly/1pADacO, and http://bit.ly/1pEpXjo
    
    # Priors for random coefficients:
    to_vector(B_g_std) ~ normal(0, 10); // implies: B_g ~ multi_normal(gamma_B_g, CovarianceMatrix);

    # Hyperpriors
    gamma_B_g ~ normal(0, 5);
    sigma_B_g_sigma ~ cauchy(0, 2.5);
    L_rho_B_g ~ lkj_corr_cholesky(3);

    //########################################################################
    // Model missing data. Need to do this explicitly in Stan.
    dbh_miss ~ normal(0, 10);

    //########################################################################
    // Main likelihood
    for (bl_num in 1:n_blocks) {
        n_rows <- block_end_row[bl_num] - block_start_row[bl_num] + 1;
        { // block to allow local dbh_pred and ID vectors of varying size
            // minus 1 below since can't make dbh_pred for first period
            matrix[n_rows, block_n_periods[bl_num] - 1] dbh_pred;
            int row_overall;
            int this_tree_ID;
            int this_plot_ID;
            int this_site_ID;
            int this_genus_ID;

            for (row_bl in 1:n_rows) {
                // convenience counter to track position in full dataset
                row_overall <- row_bl + block_start_row[bl_num] - 1;

                this_tree_ID <- tree_ID[row_overall];
                this_plot_ID <- plot_ID[this_tree_ID];
                this_site_ID <- site_ID[this_tree_ID];
                this_genus_ID <- genus_ID[this_tree_ID];

                for (t in 2:(block_n_periods[bl_num])) {
                    // minus 1 is because dbh_pred is one column smaller than 
                    // dbh_latent matrix
                    dbh_pred[row_bl, t - 1] <- B[1] * WD[this_tree_ID] +
                        B[2] * WD_sq[this_tree_ID] +
                        int_ijk[this_tree_ID] +
                        int_jk[this_plot_ID] +
                        int_k[this_site_ID] +
                        B_g[this_genus_ID, 1] +
                        B_g[this_genus_ID, 2] * precip[row_overall, t] +
                        B_g[this_genus_ID, 3] * precip_sq[row_overall, t] +
                        B_g[this_genus_ID, 4] * (B_T[this_site_ID, 1] + B_T[this_site_ID, 2] * temp[row_overall, t] + B_T[this_site_ID, 3] * elev[this_plot_ID]) +
                        B_g[this_genus_ID, 5] * square(B_T[this_site_ID, 1] + B_T[this_site_ID, 2] * temp[row_overall, t] + B_T[this_site_ID, 3] * elev[this_plot_ID]) +
                        B_g[this_genus_ID, 6] * dbh_latent[row_overall, t] +
                        B_g[this_genus_ID, 7] * dbh_latent_sq[row_overall, t];
                }
            }

            // Model dbh_latent
            to_vector(block(dbh_latent, block_start_row[bl_num], 1, block_end_row[bl_num], block_n_periods[bl_num])) ~ normal(to_vector(dbh_pred), sigma_proc);
            //print("Predicted dbh: ", dbh_pred);
        }
        // Model dbh with mean equal to latent dbh
        to_vector(block(dbh, block_start_row[bl_num], 2, block_end_row[bl_num], block_n_periods[bl_num] - 1)) ~ normal(to_vector(block(dbh_latent, block_start_row[bl_num], 2, block_end_row[bl_num], block_n_periods[bl_num] - 1)), sigma_obs);
    }

    # Specially handle first latent dbh
    to_vector(block(dbh_latent, 1, 1, n_tree, 1)) ~ normal(0, 10);
}
