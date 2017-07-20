#ifndef GUARD_utils_h
#define GUARD_utils_h

#include <math.h>

const int ITMAX = 1001; //added 1 to allow <= to be <
const double EPS = 3.0e-7;
const double FPMIN = 1.0e-30;
const double INCORRECT_INPUT = -1.0;
const double DOES_NOT_CONVERGE =  2 * INCORRECT_INPUT;
const int MAXIT=101; //for betacf; added 1 to allow <= to be <

// #include<Rcpp.h>

//functions
double timesTwo(double);
double gamma(double, double);
double gammaScaled(double, double);
double gammaln(double);
double erf(double);
double normcdf(double, double, double);
double betacf(double, double, double);
double betai(double, double, double);
double tcdf(double, double);

#endif
