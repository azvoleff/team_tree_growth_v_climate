#include<Rcpp.h>

using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::NumericMatrix left_align(Rcpp::NumericMatrix x, int indent=0) {
    NumericMatrix out(x.nrow(), x.ncol());

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
        int out_col = indent;
        for (int j=first_obs_period; j <= last_obs_period; j++) {
            out(i, out_col) = x(i, j);
            out_col++;
        }

        // Fill NAs in out for all periods after the last observation, and for 
        // all periods up until indent
        for (int j=0; j < indent; j++) {
            out(i, j) = NA_REAL;
        }
        for (int j=out_col; j < x.ncol(); j++) {
            out(i, j) = NA_REAL;
        }
    }

    return(out);
}
