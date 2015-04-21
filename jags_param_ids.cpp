#include <Rcpp.h>

#include <string>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]

const boost::regex re_scalar("^([a-zA-Z_]*)$");
const boost::regex re_vector("^([a-zA-Z_]*)\\[([0-9]*)\\]$");
const boost::regex re_matrix_2d("^([a-zA-Z_]*)\\[([0-9]*),([0-9]*)\\]$");

// [[Rcpp::export]]
DataFrame jags_param_ids(CharacterVector s) {
    int n = s.size();
    
    std::vector<std::string> param(n);
    NumericVector row(n);
    NumericVector col(n);
    
    for (int i=0; i<n; i++) {
        boost::match_results<std::string::const_iterator> matches;
        if (boost::regex_match(as<std::string>(s[i]), matches, re_scalar)) {
            param[i] = matches[1];
            row[i] = 0;
            col[i] = 0;
        } else if (boost::regex_match(as<std::string>(s[i]), matches, re_vector)) {
            param[i] = matches[1];
            row[i] = boost::lexical_cast<long>(matches[2]);
            col[i] = 0;
        } else if (boost::regex_match(as<std::string>(s[i]), matches, re_matrix_2d)) {
            param[i] = matches[1];
            row[i] = boost::lexical_cast<long>(matches[2]);
            col[i] = boost::lexical_cast<long>(matches[3]);
        } else {
            stop("could not identify parameter type for position " + boost::lexical_cast<std::string>(i + 1));
        }
    }

    return DataFrame::create(Named("Parameter_Base") = param,
                             Named("row_ID") = row,
                             Named("col_ID") = col);
}
