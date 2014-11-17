data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> n_genus;
    int<lower=0> n_B;
    int<lower=0> n_B_g;
    int<lower=0> n_B_T;
    int<lower=0> n_period;
    real<lower=0> sigma_obs_lower;
    int<lower=0> tree_ID[n_tree];
    int<lower=0> plot_ID[n_tree];
    int<lower=0> site_ID[n_tree];
    int<lower=0> genus_ID[n_tree];
    int<lower=0> n_miss;
    int<lower=0> n_obs;
    int<lower=0> obs_indices_tree[n_obs];
    int<lower=0> obs_indices_period[n_obs];
    int<lower=0> miss_indices_tree[n_miss];
    int<lower=0> miss_indices_period[n_miss];
    int<lower=0> n_blocks;
    int<lower=0> bl_size[n_blocks];
    int<lower=0> bl_st[n_blocks];
    int<lower=0> bl_end[n_blocks];
    vector[n_obs] first_obs_period;
    vector[n_obs] dbh_obs;
    vector[n_tree] WD;
    vector[n_tree] elev;
    matrix[n_tree, n_period + 1] temp;
    matrix[n_tree, n_period + 1] precip;
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
    vector[n_period * n_site] int_t_std;
    real<lower=sigma_obs_lower> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_int_ijk;
    real<lower=0> sigma_int_jk;
    real<lower=0> sigma_int_k;
    real<lower=0> sigma_int_t;
}

transformed parameters {
    matrix[n_tree, n_period + 1] dbh;
    // Don't need + 1 for dbh_pred since can't predict dbh for first period
    matrix[n_tree, n_period] dbh_pred;
    vector[n_tree] int_ijk;
    vector[n_plot] int_jk;
    vector[n_site] int_k;
    vector[n_period * n_site] int_t;
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
    int_t <- sigma_int_t * int_t_std; // int_t ~ normal(0, sigma_int_t)

    B_g <- transpose(rep_matrix(gamma_B_g, n_genus) + diag_pre_multiply(sigma_B_g_sigma, L_rho_B_g) * B_g_std);

    { // block to allow local variables
        real temp_model;
        int this_tree_ID;
        for (bl_num in 1:n_blocks) {
            for (row_num in bl_st[bl_num]:bl_end[bl_num]) {
                this_tree_ID <- tree_ID[row_num];
                for (t in 2:(bl_size[bl_num])) {
                    temp_model <- temp[row_num, t] +
                        B_T[site_ID[this_tree_ID], 1] * 
                        B_T[site_ID[this_tree_ID], 2] * elev[plot_ID[this_tree_ID]];
                    // minus 1 is because dbh_pred is one column smaller than 
                    // dbh_latent matrix
                    dbh_pred[row_num, t - 1] <- B[1] * WD[this_tree_ID] +
                        B[2] * square(WD[this_tree_ID]) +
                        int_ijk[this_tree_ID] +
                        int_jk[plot_ID[this_tree_ID]] +
                        int_k[site_ID[this_tree_ID]] +
                        int_t[(first_obs_period[this_tree_ID] + t - 2) * site_ID[this_tree_ID]] +
                        B_g[genus_ID[this_tree_ID], 1] +
                        B_g[genus_ID[this_tree_ID], 2] * precip[row_num, t] +
                        B_g[genus_ID[this_tree_ID], 3] * square(precip[row_num, t]) +
                        B_g[genus_ID[this_tree_ID], 4] * temp_model +
                        B_g[genus_ID[this_tree_ID], 5] * square(temp_model) +
                        B_g[genus_ID[this_tree_ID], 6] * dbh_latent[row_num, t - 1] +
                        B_g[genus_ID[this_tree_ID], 7] * square(dbh_latent[row_num, t - 1]);
                    // print("latent=", dbh_latent[row_num, t - 1], ", pred=", dbh_pred[row_num, t - 1]);
                }
            }
        }
    }
}

model {
    //########################################################################
    // Fixed effects
    B ~ normal(0, 10);

    //########################################################################
    // Temperature model (per site)
    to_vector(B_T) ~ normal(0, 10);

    //########################################################################
    // Observation and process error
    sigma_obs ~ cauchy(sigma_obs_lower, 1);
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
    # Period random effects (crossed)
    int_t_std ~ normal(0, 1); // Matt trick
    sigma_int_t ~ cauchy(0, 1);

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
    
    // Specially handle first latent dbh
    to_vector(block(dbh_latent, 1, 1, n_tree, 1)) ~ normal(0, 10);
    { // block to allow local variable
        int bl_rows;
        for (bl_num in 1:n_blocks) {
            bl_rows <- bl_end[bl_num] - bl_st[bl_num] + 1;
            // Model dbh_latent. Don't model first column as there is no
            // prediction for the first dbh (first column is handled above).
            to_vector(block(dbh_latent, bl_st[bl_num], 2, bl_rows, bl_size[bl_num] - 1)) ~ normal(to_vector(block(dbh_pred, bl_st[bl_num], 1, bl_rows, bl_size[bl_num] - 1)), sigma_proc);
            // Model dbh with mean equal to latent dbh.
            to_vector(block(dbh, bl_st[bl_num], 1, bl_rows, bl_size[bl_num])) ~ normal(to_vector(block(dbh_latent, bl_st[bl_num], 1, bl_rows, bl_size[bl_num])), sigma_obs);
        }
    }
}
