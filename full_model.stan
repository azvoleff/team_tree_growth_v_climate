data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> n_genus;
    int<lower=0> n_period;
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
    vector[n_obs] dbh_obs;
    vector[n_tree] WD;
    matrix[n_tree, max_obs_per_tree] spi;
}

parameters {
    matrix[n_tree, max_obs_per_tree] dbh_latent;
    vector[n_miss] dbh_miss;
    real<lower=-10, upper=10> intercept;
    real<lower=-10, upper=10> slp_dbh;
    real<lower=-10, upper=10> slp_dbh_sq;
    real<lower=-10, upper=10> slp_WD;
    real<lower=-10, upper=10> slp_WD_sq;
    real<lower=-10, upper=10> slp_spi;
    real<lower=-10, upper=10> inter_spi_WD;
    real<lower=-10, upper=10> inter_spi_dbh;
    vector[n_tree] b_ijk_std;
    vector[n_plot] b_jk_std;
    vector[n_site] b_k_std;
    vector[n_genus] b_g_std;
    vector[n_period] b_t_std;
    real<lower=0.003566014> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0> sigma_ijk;
    real<lower=0> sigma_jk;
    real<lower=0> sigma_k;
    real<lower=0> sigma_t;
    real<lower=0> sigma_g;
}

transformed parameters {
    vector[n_tree] WD_sq;
    matrix[n_tree, max_obs_per_tree] dbh;
    vector[n_tree] b_ijk;
    vector[n_plot] b_jk;
    vector[n_site] b_k;
    vector[n_genus] b_g;
    vector[n_period] b_t;

    for (n in 1:n_miss) {
        dbh[miss_indices_tree[n], miss_indices_period[n]] <- dbh_miss[n];
    }

    for (n in 1:n_obs) {
        dbh[obs_indices_tree[n], obs_indices_period[n]] <- dbh_obs[n];
    }

    for (n in 1:n_tree) {
        WD_sq[n] <- pow(WD[n], 2);
    }

    // Matt trick- see http://bit.ly/1qz4NC6
    b_ijk <- sigma_ijk * b_ijk_std; // b_ijk ~ normal(0, sigma_ijk)
    b_jk <- sigma_jk * b_jk_std; // b_jk ~ normal(0, sigma_jk)
    b_k <- sigma_k * b_k_std; // b_k ~ normal(0, sigma_k)
    b_g <- sigma_g * b_g_std; // b_g ~ normal(0, sigma_g)
    b_t <- sigma_t * b_t_std; // b_t ~ normal(0, sigma_t)
}

model {
    dbh_miss ~ normal(0, 10);

    for (i in 1:n_tree) {
        for (t in (first_obs_period[i]):last_obs_period[i]) {
            dbh[i, t] ~ normal(dbh_latent[i, t], sigma_obs);
        }

        # Specially handle first latent dbh (following Eitzen, 2013):
        dbh_latent[i, first_obs_period[i]] ~ normal(0, 10);

        for (t in (first_obs_period[i] + 1):last_obs_period[i]) {
            {
                real dbh_predicted;
                //print(i, ", ", t);
                dbh_predicted <- intercept +
                    slp_dbh * dbh_latent[i, t - 1] +
                    slp_dbh_sq * pow(dbh_latent[i, t - 1], 2) +
                    slp_WD * WD[i] +
                    slp_WD_sq * WD_sq[i] +
                    slp_spi * spi[i, t] +
                    inter_spi_dbh * spi[i, t] * dbh_latent[i, t - 1] +
                    inter_spi_WD * spi[i, t] * WD[i] +
                    b_ijk[i] +
                    b_jk[plot_ID[i]] +
                    b_k[site_ID[i]] +
                    b_t[t] +
                    b_g[genus_ID[i]];

                dbh_latent[i, t] ~ normal(dbh_predicted, sigma_proc);

                /* print("Latent dbh (t): ", dbh_latent[i, t]); */
                /* print("Latent dbh (t - 1): ", dbh_latent[i, t - 1]); */
                /* print("Predicted dbh (t): ", dbh_predicted); */
            }
        }
    }

    b_ijk_std ~ normal(0, 1); // Matt trick
    b_jk_std ~ normal(0, 1); // Matt trick
    b_k_std ~ normal(0, 1); // Matt trick
    b_g_std ~ normal(0, 1); // Matt trick
    b_t_std ~ normal(0, 1); // Matt trick

    sigma_obs ~ cauchy(0, 10);
    sigma_proc ~ cauchy(0, 10);
    sigma_ijk ~ cauchy(0, 10);
    sigma_jk ~ cauchy(0, 10);
    sigma_k ~ cauchy(0, 10);
    sigma_t ~ cauchy(0, 10);
    sigma_g ~ cauchy(0, 10);
}
