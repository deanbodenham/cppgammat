#ifndef GUARD_utilsRcpp_h
#define GUARD_utilsRcpp_h

#include "utils.h"
#include<Rcpp.h>

Rcpp::NumericVector timesTwoRcpp(Rcpp::NumericVector);

Rcpp::NumericVector gammaRcpp(Rcpp::NumericVector);

Rcpp::NumericVector gammaScaledRcpp(Rcpp::NumericVector);

Rcpp::NumericVector gammalnRcpp(Rcpp::NumericVector);

Rcpp::NumericVector erfRcpp(Rcpp::NumericVector);

Rcpp::NumericVector normcdfRcpp(Rcpp::NumericVector, 
                                Rcpp::NumericVector, 
                                Rcpp::NumericVector);

Rcpp::NumericVector betacfRcpp(Rcpp::NumericVector, 
                               Rcpp::NumericVector, 
                               Rcpp::NumericVector);

Rcpp::NumericVector betaiRcpp(Rcpp::NumericVector, 
                              Rcpp::NumericVector, 
                              Rcpp::NumericVector);

Rcpp::NumericVector tcdfRcpp(Rcpp::NumericVector, 
                             Rcpp::NumericVector);

#endif
