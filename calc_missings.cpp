#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List calc_missings(Rcpp::NumericMatrix x) {
    // Make two matrices that will represent missing and observed values in x 
    // in long format (where the first row is the row number in the original 
    // matrix x, and the second column is the period number column number in 
    // the original matrix x)
    //
    // Note: allocating the two arrays the way I do in the next two lines is 
    // lazy.  Could optimize this to not need to allocate so much memory up 
    // front by dynamically adding rows in blocks.
    NumericMatrix miss(x.nrow() * x.ncol(), 2);
    NumericMatrix obs(x.nrow() * x.ncol(), 2);

    // Counters to track how many rows of the above matrices are used
    int num_miss = 0;
    int num_obs = 0;
    
    for (int i=0; i < x.nrow(); i++) {
        // First calculate the first and last observation period (as missing 
        // data only counts if it is within the period of the observed data)
        int first_obs_period = x.ncol();
        int last_obs_period = 0;
        for (int j=0; j < x.ncol(); j++) {
            if (j < first_obs_period && !std::isnan(x(i, j))) first_obs_period = j;
            if (j > last_obs_period && !std::isnan(x(i, j))) last_obs_period = j;
        }

        // Now calculate the observed and missing data points within the bounds 
        // of the observed data
        for (int j=first_obs_period; j <= last_obs_period; j++) {
            if (std::isnan(x(i, j))) {
                miss(num_miss, 0) = i + 1; // R indices start a 1
                miss(num_miss, 1) = j + 1;
                num_miss++;
            } else {
                obs(num_obs, 0) = i + 1; // R indices start a 1
                obs(num_obs, 1) = j + 1;
                num_obs++;
            }
        }
    }

    // Don't return empty, unused rows in obs and miss. Return NULL for obs or 
    // miss if either is empty.
    List out;
    if (num_miss == 0) {
        out["miss"] = R_NilValue;
    } else {
        NumericMatrix miss_out(num_miss, miss.ncol());
        for (int i=0; i < num_miss; i++) {
            miss_out.row(i) = miss.row(i);
        }
        out["miss"] = miss_out;
    }
    if (num_obs == 0) {
        out["obs"] = R_NilValue;
    } else {
        NumericMatrix obs_out(num_obs, obs.ncol());
        for (int i=0; i < num_obs; i++) {
            obs_out.row(i) = obs.row(i);
        }
        out["obs"] = obs_out;
    }

    return(out);
}
