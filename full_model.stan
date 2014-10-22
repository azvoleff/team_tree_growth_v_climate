data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> n_genus;
    int<lower=0> n_period;
    int<lower=0> n_B;
    int<lower=0> n_B_g;
    int<lower=0> plot_ID[n_tree];
    int<lower=0> site_ID[n_tree];
    int<lower=0> genus_ID[n_tree];
    int<lower=0> max_obs_per_tree;
    int<lower=0> first_obs_period[n_tree];
    int<lower=0> last_obs_period[n_tree];
    int<lower=0> n_miss;
    int<lower=0> n_obs;
    int<lower=0> obs_indices_tree[n_obs];
    int<lower=0> obs_indices_period[n_obs];
    int<lower=0> miss_indices_tree[n_miss];
    int<lower=0> miss_indices_period[n_miss];
    matrix[n_B_g, n_B_g] W;
    vector[n_obs] dbh_obs;
    vector[n_tree] WD;
    vector[n_tree] WD_sq;
    matrix[n_tree, max_obs_per_tree] mcwd;
    matrix[n_tree, max_obs_per_tree] mcwd_sq;
}

parameters {
    matrix[n_tree, max_obs_per_tree] dbh_latent; 
    vector[n_B_g] B_g_raw_temp;
    vector[n_miss] dbh_miss;
    vector[n_B] B;
    vector[n_B_g] mu_B_g_raw;
    vector[n_B_g] xi;
    corr_matrix[n_B_g] Tau_B_g_raw;
    vector[n_tree] int_ijk_std;
    vector[n_plot] int_jk_std;
    vector[n_site] int_k_std;
    vector[n_period] int_t_std;
    real<lower=0.003566014> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_int_ijk;
    real<lower=0> sigma_int_jk;
    real<lower=0> sigma_int_k;
    real<lower=0> sigma_int_t;
}

transformed parameters {
    matrix[n_tree, max_obs_per_tree] dbh;
    vector[n_tree] int_ijk;
    vector[n_plot] int_jk;
    vector[n_site] int_k;
    vector[n_period] int_t;
    real<lower=0> df;

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

    df <- n_B_g + 1;
}

model {
    matrix[n_tree, max_obs_per_tree] dbh_predicted;
    matrix[n_genus, n_B_g] B_g;
    matrix[n_genus, n_B_g] B_g_raw;
    matrix[n_B_g, n_B_g] rho_B_g;
    matrix[n_B_g, n_B_g] Sigma_B_g_raw;
    vector[n_B_g] mu_B_g;
    vector[n_B_g] sigma_B_g;

    ##########################################################################
    # Fixed effects
    B ~ normal(0, 10);

    ##########################################################################
    # Observation and process error
    sigma_obs ~ cauchy(0, 10);
    sigma_proc ~ cauchy(0, 10);

    ##########################################################################
    # Nested random effects
    int_ijk_std ~ normal(0, 1); // Matt trick
    sigma_int_ijk ~ cauchy(0, 10);

    int_jk_std ~ normal(0, 1); // Matt trick
    sigma_int_jk ~ cauchy(0, 10);

    sigma_int_k ~ cauchy(0, 10);
    int_k_std ~ normal(0, 1); // Matt trick

    ##########################################################################
    # Period random effects (crossed)
    int_t_std ~ normal(0, 1); // Matt trick
    sigma_int_t ~ cauchy(0, 10);

    ##########################################################################
    # Correlated random effects at genus level (crossed). Based on Gelman and 
    # Hill pg. 377-378, and http://bit.ly/1pADacO
    mu_B_g_raw ~ normal(0, 100);
    xi ~ uniform(0, 100);

    mu_B_g <- xi .* mu_B_g_raw;

    Tau_B_g_raw ~ wishart(df, W);
    Sigma_B_g_raw <- inverse(Tau_B_g_raw);

    for (j in 1:n_genus) {
        B_g_raw_temp ~ multi_normal(mu_B_g_raw, Sigma_B_g_raw);
        for (k in 1:n_B_g) {
            B_g_raw[j, k] <- B_g_raw_temp[k];
            B_g[j, k] <- xi[k] * B_g_raw[j, k];
        }
    }

    for (k in 1:n_B_g) {
        for (k_prime in 1:n_B_g) {
            rho_B_g[k, k_prime] <- Sigma_B_g_raw[k, k_prime] / sqrt(Sigma_B_g_raw[k, k] * Sigma_B_g_raw[k_prime, k_prime]);
        }
        sigma_B_g[k] <- fabs(xi[k]) * sqrt(Sigma_B_g_raw[k,k]);
    }

    ##########################################################################
    # Model missing data. Need to do this explicitly in Stan
    dbh_miss ~ normal(0, 10);

    ##########################################################################
    # Main likelihood
    for (i in 1:n_tree) {
        //print("Starting tree ", i);
        # Specially handle first latent dbh (following Eitzen, 2013):
        dbh_latent[i, first_obs_period[i]] ~ normal(0, 10);

        for (t in (first_obs_period[i] + 1):last_obs_period[i]) {
            //print(i, ", ", t);
            dbh_predicted[i, t] <- B[1] +
                B[2] * mcwd[i, t - 1] +
                B[3] * mcwd_sq[i, t - 1] +
                B[4] * dbh_latent[i, t - 1] +
                B[5] * pow(dbh_latent[i, t - 1], 2) +
                B[6] * WD[i] +
                B[7] * WD_sq[i] +
                int_ijk[i] +
                int_jk[plot_ID[i]] +
                int_k[site_ID[i]] +
                int_t[t - 1] +
                B_g[genus_ID[i], 1] +
                B_g[genus_ID[i], 2] * mcwd[i, t - 1] +
                B_g[genus_ID[i], 3] * mcwd_sq[i, t - 1] +
                B_g[genus_ID[i], 4] * dbh_latent[i, t - 1] +
                B_g[genus_ID[i], 5] * pow(dbh_latent[i, t - 1], 2);
            dbh_latent[i, t] ~ normal(dbh_predicted[i, t], sigma_proc);
            /* print("Latent dbh (t): ", dbh_latent[i, t]); */
            /* print("Latent dbh (t - 1): ", dbh_latent[i, t - 1]); */
            /* print("Predicted dbh (t): ", dbh_predicted); */
        }

        for (t in (first_obs_period[i]):last_obs_period[i]) {
            dbh[i, t] ~ normal(dbh_latent[i, t], sigma_obs);
        }
    }
}
