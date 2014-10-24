data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> n_genus;
    int<lower=0> n_B;
    int<lower=0> n_B_g;
    int<lower=0> n_period;
    int<lower=0> plot_ID[n_tree];
    int<lower=0> site_ID[n_tree];
    int<lower=0> genus_ID[n_tree];
    int<lower=0> first_obs_period[n_tree];
    int<lower=0> last_obs_period[n_tree];
    int<lower=0> n_miss;
    int<lower=0> n_obs;
    int<lower=0> obs_indices_tree[n_obs];
    int<lower=0> obs_indices_period[n_obs];
    int<lower=0> miss_indices_tree[n_miss];
    int<lower=0> miss_indices_period[n_miss];
    vector[n_obs] dbh_obs;
    vector[n_tree] WD;
    vector[n_tree] WD_sq;
    matrix[n_tree, n_period + 1] mcwd;
    matrix[n_tree, n_period + 1] mcwd_sq;
}

parameters {
    matrix[n_tree, n_period + 1] dbh_latent; 
    vector[n_miss] dbh_miss;
    vector[n_B] B;
    matrix[n_B_g, n_genus] B_g_std;
    vector[n_B_g] gamma_B_g;
    cholesky_factor_corr[n_B_g] L_rho_B_g;
    vector<lower=0>[n_B_g] sigma_B_g_sigma;
    vector[n_tree] int_ijk_std;
    vector[n_plot] int_jk_std;
    vector[n_site] int_k_std;
    //vector[n_period] int_t_std;
    real<lower=0.0002537298> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_int_ijk;
    real<lower=0> sigma_int_jk;
    real<lower=0> sigma_int_k;
    //real<lower=0> sigma_int_t;
}

transformed parameters {
    matrix[n_tree, n_period + 1] dbh;
    vector[n_tree] int_ijk;
    vector[n_plot] int_jk;
    vector[n_site] int_k;
    //vector[n_period] int_t;
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
    //int_t <- sigma_int_t * int_t_std; // int_t ~ normal(0, sigma_int_t)

    B_g <- transpose(rep_matrix(gamma_B_g, n_genus) + diag_pre_multiply(sigma_B_g_sigma, L_rho_B_g) * B_g_std);
}

model {
    //matrix[n_tree, n_period + 1] dbh_predicted;
    real dbh_predicted;

    //########################################################################
    // Fixed effects
    B ~ normal(0, 10);

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
    # Period random effects (crossed)
    //int_t_std ~ normal(0, 1); // Matt trick
    //sigma_int_t ~ cauchy(0, 1);

    //########################################################################
    // Correlated random effects at genus level (crossed). Based on Gelman and 
    // Hill pg. 377-378, and http://bit.ly/1pADacO, and http://bit.ly/1pEpXjo
    
    # Priors for random coefficients:
    to_vector(B_g_std) ~ normal(0, 1); // implies: B_g ~ multi_normal(gamma_B_g, CovarianceMatrix);

    # Hyperpriors
    gamma_B_g ~ normal(0, 5);
    sigma_B_g_sigma ~ cauchy(0, 2.5);
    L_rho_B_g ~ lkj_corr_cholesky(3);

    //########################################################################
    // Model missing data. Need to do this explicitly in Stan
    dbh_miss ~ normal(0, 10);

    //########################################################################
    // Main likelihood
    for (i in 1:n_tree) {
        //print("Starting tree ", i);
        # Specially handle first latent dbh (following Eitzen, 2013):
        dbh_latent[i, first_obs_period[i]] ~ normal(0, 10);

        for (t in (first_obs_period[i] + 1):last_obs_period[i]) {
            //print(i, ", ", t);
            dbh_predicted <- B[1] +
                B[2] * mcwd[i, t - 1] +
                B[3] * mcwd_sq[i, t - 1] +
                B[4] * dbh_latent[i, t - 1] +
                B[5] * pow(dbh_latent[i, t - 1], 2) +
                B[6] * WD[i] +
                B[7] * WD_sq[i] +
                int_ijk[i] +
                int_jk[plot_ID[i]] +
                int_k[site_ID[i]] +
                B_g[genus_ID[i], 1] +
                B_g[genus_ID[i], 2] * mcwd[i, t - 1] +
                B_g[genus_ID[i], 3] * mcwd_sq[i, t - 1] +
                B_g[genus_ID[i], 4] * dbh_latent[i, t - 1] +
                B_g[genus_ID[i], 5] * pow(dbh_latent[i, t - 1], 2);
                //int_t[t - 1] +
            dbh_latent[i, t] ~ normal(dbh_predicted, sigma_proc);
            /* print("Latent dbh (t): ", dbh_latent[i, t]); */
            /* print("Latent dbh (t - 1): ", dbh_latent[i, t - 1]); */
            /* print("Predicted dbh (t): ", dbh_predicted); */
        }

        for (t in (first_obs_period[i]):last_obs_period[i]) {
            dbh[i, t] ~ normal(dbh_latent[i, t], sigma_obs);
        }
    }
}
