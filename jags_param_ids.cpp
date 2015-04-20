#include <Rcpp.h>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]

#include <string>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>

bool validate_vector(const std::string& s) {
   static const boost::regex e("\\A[a-zA-Z_]*\\[[0-9]*\\]");
   return boost::regex_match(s, e);
}

bool validate_matrix_2d(const std::string& s) {
   static const boost::regex e("\\A[a-zA-Z_]*\\[[0-9]*,[0-9]*\\]");
   return boost::regex_match(s, e);
}

const boost::regex re_param_name("\\A([a-zA-Z_]*)");
const boost::regex re_vector("\\A([a-zA-Z_]*)\\[([0-9]*)\\]");
const boost::regex re_matrix_2d("\\A([a-zA-Z_]*)\\[([0-9]*),([0-9]*)\\]");

const std::string param_name("\\1");
const std::string vector_pos("\\2");
const std::string matrix_2d_row("\\2");
const std::string matrix_2d_col("\\3");

std::string get_param_name(const std::string& s) {
   return boost::regex_replace(s, re_param_name, "\\1", boost::match_default | boost::format_sed);
}

std::string get_vector_pos(const std::string& s) {
   return boost::regex_replace(s, re_vector, vector_pos, boost::match_default | boost::format_sed);
}

std::string get_matrix_2d_row(const std::string& s) {
   return boost::regex_replace(s, re_matrix_2d, matrix_2d_row, boost::match_default | boost::format_sed);
}

std::string get_matrix_2d_col(const std::string& s) {
   return boost::regex_replace(s, re_matrix_2d, matrix_2d_col, boost::match_default | boost::format_sed);
}

// [[Rcpp::export]]
Rcpp::DataFrame jags_param_ids(std::vector<std::string> s) {
    int n = s.size();
    
    std::vector<std::string> param(n);
    Rcpp::NumericVector row(n);
    Rcpp::NumericVector col(n);
    
    for (int i=0; i<n; i++) {
        param[i] = get_param_name(s[i]);

        if (validate_vector(s[i])) {
            row[i] = boost::lexical_cast<short>(get_vector_pos(s[i]));
            col[i] = NA_REAL;
        } else if (validate_matrix_2d(s[i])) {
            row[i] = boost::lexical_cast<short>(get_matrix_2d_row(s[i]));
            col[i] = boost::lexical_cast<short>(get_matrix_2d_col(s[i]));
        } else {
            row[i] = NA_REAL;
            col[i] = NA_REAL;
        }
    }
    return Rcpp::DataFrame::create(Rcpp::Named("param") = param,
                                   Rcpp::Named("row") = row,
                                   Rcpp::Named("col") = col);
}
