data {
    int<lower=0> n_tree;
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
}

parameters {
    matrix<lower=7.5, upper=300>[n_tree, max_obs_per_tree] dbh_latent;
    vector<lower=7.5, upper=300>[n_miss] dbh_miss;
    real<lower=-100, upper=100> inter;
    real<lower=-100, upper=100> slp_dbh;
    real<lower=-100, upper=100> slp_dbh_sq;
    real<lower=0> sigma_obs;
    real<lower=0> sigma_proc;
}

transformed parameters {
    matrix[n_tree, max_obs_per_tree] dbh;
    /* for (n in 1:n_miss) { */
    /*     dbh[miss_indices_tree[n], miss_indices_period[n]] ~ normal(miss_indices_tree[n], miss_indices_period[n] - 1, 10); */
    /* } */
    for (n in 1:n_miss) {
        dbh[miss_indices_tree[n], miss_indices_period[n]] <- dbh_miss[n];
    }
    for (n in 1:n_obs) {
        dbh[obs_indices_tree[n], obs_indices_period[n]] <- dbh_obs[n];
    }
}

model {
    for (tree_num in 1:n_tree) {
        # Specially handle first latent dbh (following Eitzen, 2013):
        dbh_latent[tree_num, first_obs_period[tree_num]] ~ normal(21, 100);
        # Make sure first dbh is modeled:
        dbh[tree_num, first_obs_period[tree_num]] ~ normal(dbh_latent[tree_num, first_obs_period[tree_num]], sigma_obs);
        for (obs_num in (first_obs_period[tree_num] + 1):last_obs_period[tree_num]) {
            {
                real dbh_predicted;
                dbh_predicted <- inter + slp_dbh * dbh_latent[tree_num, obs_num - 1] + slp_dbh_sq * pow(dbh_latent[tree_num, obs_num - 1], 2);
                dbh[tree_num, obs_num] ~ normal(dbh_latent[tree_num, obs_num], sigma_obs);
                dbh_latent[tree_num, obs_num] ~ normal(dbh_predicted, sigma_proc);
                /* print("dbh_predicted=", dbh_predicted[tree_num, obs_num], */
                /*       ", dbh_latent=", dbh_latent[tree_num, obs_num], */
                /*       ", dbh=", dbh[tree_num, obs_num]) */
            }
        }
    }
    sigma_obs ~ cauchy(0, 10);
    sigma_proc ~ cauchy(0, 10);
}
