#ifndef GUARD_utilsRcpp_cpp
#define GUARD_utilsRcpp_cpp

#include "utilsRcpp.h"
#include<Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector timesTwoRcpp(Rcpp::NumericVector x_){
    //cast Rcpp NumericVector to double
    double x = Rcpp::as<double> (x_);

    //compute function
    double ans = timesTwo(x);

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector gammalnRcpp(Rcpp::NumericVector x_){
    double x = Rcpp::as<double> (x_);

    //compute function
    double ans = gammaln(x);

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}



// [[Rcpp::export]]
Rcpp::NumericVector gammaRcpp(Rcpp::NumericVector a_, Rcpp::NumericVector x_){
    //cast Rcpp NumericVectors to doubles
    double a = Rcpp::as<double> (a_);
    double x = Rcpp::as<double> (x_);

    double ans;
    //compute function
    //first check bad case
    if (x < 0.0 || a <= 0.0){
        //error, stop with error message
        Rcpp::stop("Both x and a must be positive.");

        //this return is not reached
        return NA_REAL;             
    } 


    //if it has made it this far, then all is good
    ans = gamma(a, x);
    if (ans < 0.0){

        //there is an error, so check the error
        //weird scaling to check in between INCORRECT_INPUT 
        //and DOES_NOT_CONVERGE
        if (ans < (INCORRECT_INPUT * 1.5)){
            //now must be DOES_NOT_CONVERGE
            Rcpp::stop("WARNING: Does not converge");
        } else {
            Rcpp::stop("ERROR: incorrect input - both x and a must be positive.");
        }
    
    }

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}



// [[Rcpp::export]]
Rcpp::NumericVector gammaScaledRcpp(Rcpp::NumericVector x_, Rcpp::NumericVector a_){
    //SCALED VERSION, same as R's pgamma (including same order of arguments)
    //cast Rcpp NumericVectors to doubles
    double a = Rcpp::as<double> (a_);
    double x = Rcpp::as<double> (x_);

    double ans;
    //compute function
    //first check bad case
    if (x < 0.0 || a <= 0.0){
        //error, stop with error message
        Rcpp::stop("Both x and a must be positive.");

        //this return is not reached
        return NA_REAL;             
    } 


    //if it has made it this far, then all is good
    ans = gammaScaled(x, a);
    if (ans < 0.0){

        //there is an error, so check the error
        //weird scaling to check in between INCORRECT_INPUT 
        //and DOES_NOT_CONVERGE
        if (ans < (INCORRECT_INPUT * 1.5)){
            //now must be DOES_NOT_CONVERGE
            Rcpp::stop("WARNING: Does not converge");
        } else {
            Rcpp::stop("ERROR: incorrect input - both x and a must be positive.");
        }
    
    }

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector erfRcpp(Rcpp::NumericVector x_){
    //cast Rcpp NumericVector to double
    double x = Rcpp::as<double> (x_);

//     if (x < 0.0){
//         //error, stop with error message
//         Rcpp::stop("x must be strictly positive.");
//     }

    //compute function
    double ans = erf(x);

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}




// [[Rcpp::export]]
Rcpp::NumericVector normcdfRcpp(Rcpp::NumericVector x_, 
                                Rcpp::NumericVector mu_,
                                Rcpp::NumericVector sigma_){

    //cast Rcpp NumericVectors to doubles
    double x = Rcpp::as<double> (x_);
    double mu = Rcpp::as<double> (mu_);
    double sigma = Rcpp::as<double> (sigma_);

    if (sigma <= 0.0){
        //error, stop with error message
        Rcpp::stop("sigma must be strictly positive.");
    }

    //compute function
    double ans = normcdf(x, mu, sigma);

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector betacfRcpp(Rcpp::NumericVector a_,
                               Rcpp::NumericVector b_,
                               Rcpp::NumericVector x_){

    //cast Rcpp NumericVectors to doubles
    double a = Rcpp::as<double> (a_);
    double b = Rcpp::as<double> (b_);
    double x = Rcpp::as<double> (x_);


    //if it has made it this far, then all is good
    double ans = betacf(a, b, x);
    if (ans < 0.0){
        //now must be DOES_NOT_CONVERGE
        Rcpp::stop("WARNING: Does not converge");
    }

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector betaiRcpp(Rcpp::NumericVector a_,
                              Rcpp::NumericVector b_,
                              Rcpp::NumericVector x_){

    //cast Rcpp NumericVectors to doubles
    double a = Rcpp::as<double> (a_);
    double b = Rcpp::as<double> (b_);
    double x = Rcpp::as<double> (x_);

    if (x < 0.0 || x > 1.0 || a <= 0 || b <= 0){
        Rcpp::stop("ERROR: incorrect input - x must be between 0 and 1, and a and b must both be strictly positive.");
    }

    //if it has made it this far, then all is good
    double ans = betai(a, b, x);
    if (ans < 0.0){
        //now must be INCORRECT_INPUT
        Rcpp::stop("ERROR: incorrect input - x must be between 0 and 1, and a and b must both be strictly positive.");
    }

    //cast double back to Rcpp::NumericVector
    return Rcpp::NumericVector::create(ans);
}


// [[Rcpp::export]]
Rcpp::NumericVector tcdfRcpp(Rcpp::NumericVector t_, 
                             Rcpp::NumericVector nu_){

    //cast back to doubles
    double t = Rcpp::as<double> (t_);
    double nu = Rcpp::as<double> (nu_);

    if (nu <= 0){
        Rcpp::stop("ERROR: incorrect input - nu must be  strictly greater than 0.");
    }

    //check answer
    double ans = tcdf(t, nu);

    //could be incorrect input or not converging
    if (ans < 0.0){
        if (ans < (INCORRECT_INPUT * 1.5)){
            //now must be DOES_NOT_CONVERGE
            Rcpp::stop("WARNING: t cdf Does not converge");
        } else {
            //now must be INCORRECT_INPUT
            Rcpp::stop("ERROR: incorrect input for t cdf.");
        }
    }

    return Rcpp::NumericVector::create(ans);
}


#endif
