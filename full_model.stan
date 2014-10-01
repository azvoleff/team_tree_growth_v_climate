data {
    int<lower=0> n_tree;
    int<lower=0> n_plot;
    int<lower=0> n_site;
    int<lower=0> max_obs_per_tree;
    int<lower=0> first_obs_period[n_tree];
    int<lower=0> last_obs_period[n_tree];
    int<lower=0> n_miss;
    int<lower=0> n_obs;
    int<lower=0> obs_indices_tree[n_obs];
    int<lower=0> obs_indices_period[n_obs];
    int<lower=0> miss_indices_tree[n_miss];
    int<lower=0> miss_indices_period[n_miss];
    vector<lower=7.5, upper=300>[n_obs] dbh_obs;
    vector[n_tree] WD;
    matrix[n_tree, max_obs_per_tree] spi;
}

parameters {
    matrix<lower=7.5, upper=300>[n_tree, max_obs_per_tree] dbh_latent;
    vector<lower=7.5, upper=300>[n_miss] dbh_miss;
    real<lower=-100, upper=100> intercept;
    real<lower=-100, upper=100> slp_dbh;
    real<lower=-100, upper=100> slp_dbh_sq;
    real<lower=-100, upper=100> slp_WD;
    real<lower=-100, upper=100> slp_WD_sq;
    real<lower=-100, upper=100> slp_spi;
    real<lower=-100, upper=100> inter_spi_WD;
    real<lower=-100, upper=100> inter_spi_dbh;
    real<lower=0> sigma_obs;
    real<lower=0> sigma_proc;
    real<lower=0, upper=100> sigma_ijk;
    real<lower=0, upper=100> sigma_jk;
    real<lower=0, upper=100> sigma_k;
}

transformed parameters {
    vector[n_tree] WD_sq;
    matrix[n_tree, max_obs_per_tree] dbh;
    vector[n_tree] b_ijk;
    vector[n_plot] b_jk;
    vector[n_site] b_k;

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
}

model {
    for (tree_num in 1:n_tree) {
        for (obs_num in (first_obs_period[tree_num] + 1):last_obs_period[tree_num]) {
            {
                real dbh_predicted;
                dbh_predicted <- intercept +
                    slp_dbh * dbh_latent[tree_num, obs_num - 1] +
                    slp_dbh_sq * pow(dbh_latent[tree_num, obs_num - 1], 2) +
                    slp_WD * WD[tree_num] +
                    slp_WD_sq * WD_sq[tree_num] +
                    slp_spi * spi[tree_num, obs_num] +
                    inter_spi_dbh * spi[tree_num, obs_num] * dbh_latent[tree_num, obs_num - 1] +
                    inter_spi_WD * spi[tree_num, obs_num] * WD[tree_num] +
                    b_ijk[tree_num] + b_jk[plot_ID[tree_num]] + b_k[site_ID[tree_num]];

                dbh_latent[tree_num, obs_num] ~ normal(dbh_predicted, sigma_proc);
                /* print("dbh_predicted=", dbh_predicted, */
                /*       ", dbh_latent=", dbh_latent[tree_num, obs_num], */
                /*       ", dbh=", dbh[tree_num, obs_num]); */
            }
        }

        for (obs_num in (first_obs_period[tree_num]):last_obs_period[tree_num]) {
            dbh[tree_num, obs_num] ~ normal(dbh_latent[tree_num, obs_num], sigma_obs);
        }

        # Specially handle first latent dbh (following Eitzen, 2013):
        dbh_latent[tree_num, first_obs_period[tree_num]] ~ normal(0, 100);
    }
    sigma_obs ~ cauchy(0, 10);
    sigma_proc ~ cauchy(0, 10);
    b_ijk_std ~ normal(0, 1); // Matt trick
    b_jk_std ~ normal(0, 1); // Matt trick
    b_k_std ~ normal(0, 1); // Matt trick
}
