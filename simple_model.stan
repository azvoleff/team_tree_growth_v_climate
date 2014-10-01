data {
    int<lower=0> n_tree;
    int<lower=0> max_obs_per_tree;
    int<lower=0> first_obs_period[n_tree];
    int<lower=0> last_obs_period[n_tree];
    matrix[n_tree, max_obs_per_tree] log_dbh;
}

parameters {
    matrix<lower=2, upper=5.6>[n_tree, max_obs_per_tree] log_dbh_latent;
    real<lower=-3, upper=3> inter;
    real<lower=-3, upper=3> slp_dbh;
    real<lower=0> sigma_obs;
    real<lower=0> sigma_proc;
}

model {
    matrix[n_tree, max_obs_per_tree] log_growth_hat;
    matrix<lower=0>[n_tree, max_obs_per_tree] log_dbh_latent_chg;
    matrix<lower=0>[n_tree, max_obs_per_tree] growth;
    for (tree_num in 1:n_tree) {
        // Make sure first latent dbh is also modeled
        log_dbh[tree_num, first_obs_period[tree_num]] ~ normal(log_dbh_latent[tree_num, first_obs_period[tree_num]], sigma_obs);
        for (obs_num in (first_obs_period[tree_num] + 1):last_obs_period[tree_num]) {
            // Distribution of dbh, with  mean latent dbh
            log_dbh_latent_chg[tree_num, obs_num] <- (log_dbh_latent[tree_num, obs_num] - log_dbh_latent[tree_num, obs_num - 1]) / log_dbh_latent[tree_num, obs_num - 1];
            log_dbh[tree_num, obs_num] ~ normal(log_dbh_latent[tree_num, obs_num], sigma_obs);

            // Growth model
            growth[tree_num, obs_num] <- exp(log_dbh_latent[tree_num, obs_num]) - exp(log_dbh_latent[tree_num, obs_num - 1]);
            log_growth_hat[tree_num, obs_num] <- inter + slp_dbh * exp(log_dbh_latent[tree_num, obs_num - 1]);
            log(growth[tree_num, obs_num]) ~ normal(log_growth_hat[tree_num, obs_num], sigma_proc);

            // Jacobian adjustment to account for log transform of latent growth 
            // variables (log absolute determinant of transform):
            increment_log_prob(-log(fabs(growth[tree_num, obs_num])));
        }
    }

    sigma_obs ~ cauchy(0, 1);
    sigma_proc ~ cauchy(0, 1);
}
