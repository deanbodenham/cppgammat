#ifndef GUARD_utils_cpp
#define GUARD_utils_cpp

#include "utils.h"


//just multiplying a number by two
double timesTwo(double x){
    return x*2;
}


/*-------------------------------------------------------------------*/
//Gamma functions
//computation of gamma function from 
//Numerical Recipes, Chapter 6.2
/*-------------------------------------------------------------------*/


//log gamma function
//Reference: Numerical Recipes in C, Section 6.1
//Based on the Lanczos approximation
double gammaln(double xx){
    double x, y, tmp, ser;

    const double cof[6] = {76.18009172947146,
                           -86.50532032941677,
                            24.01409824083091,
                           -1.231739572450155,
                            0.1208650973866179e-2,
                           -0.5395239384953e-5};

    y=x=xx;
    tmp = x+5.5;
    tmp -= (x+0.5) * log(tmp);
    ser=1.000000000190015;
    for (int j=0; j <=5; j++){
        ser += cof[j]/++y;
    }

    return -tmp + log(2.5066282746310005*ser/x);
}




//helper function
//returns incomplete gamma function evaluated by its 
//        series representation epexpression
//returns log of Gamma(a) as gln (extra)
//
//Although incorrect input should never be reached, because it
//is already checked in gammp and gammq, we still have basic checks here
void gser(double *gamser, double a, double x, double *gln, 
          bool *doesNotConverge, bool *incorrectInput){

    //Errors are handled in Rcpp wrapper in utilsRcpp.cpp
    double sum, del, ap;
    *doesNotConverge = true;
    *incorrectInput = false;

    if (a < 0.0){
        *incorrectInput = true;
        return;
    }

    //compute the log gamma - needed later anyway
    *gln = gammaln(a);

    if (x <= 0.0){
        if (x < 0.0){
            *incorrectInput = true;
            return;
        }
        //series is zero - end of function
        *gamser=0.0;
        *doesNotConverge = false;
        return;
    } else {
        ap = a;
        del = sum = 1.0/a;
        for (int n=1; n < ITMAX; n++){
            ++ap;
            del *= x/ap;
            sum += del;
            //use fabs instead of std::abs,
            //because not forcing C++11
            if (fabs(del) < fabs(sum)*EPS){
                *gamser =  sum * exp(-x+a*log(x) - (*gln));
                *doesNotConverge = false;
                return;
            }
        } //end of for
        return;
    } //endof if else
}



//gcf
//also returns incomplete gamma function...
void gcf(double *gammcf, double a, double x, double *gln, 
         bool *doesNotConverge){
    double an, b, c, d, del, h;
    int i;

    *gln = gammaln(a);
    b = x + 1.0 - a;
    c = 1.0/FPMIN;
    d = 1.0/b;
    h = d;
    for (i=1; i < ITMAX; i++){
        an = -i * (i-a);
        b += 2.0;

        d = an*d+b;
        if (fabs(d) < FPMIN)
            d = FPMIN;

        c = b+an/c;
        if (fabs(c) < FPMIN)
            c = FPMIN;

        d = 1.0/d;
        del = d*c;
        h *= del;

        if (fabs(del-1.0) < EPS)
            break;
    }
    if (i > (ITMAX-1)){
        *doesNotConverge = true;
    }
    *gammcf = exp( -x+a*log(x) - (*gln) )*h;
    *doesNotConverge = false;
}




/*-------------------------------------------------------------------*/
//gammp function - gammq below
//returns negative values for errors
double gammp(double a, double x){
    double gamser; 
    double gammcf;
    double gln;
    bool doesNotConverge;
    bool incorrectInput;

    if ((x < 0.0) || (a <= 0.0)){
        //error - but this should never be reached because of wrapper
        return INCORRECT_INPUT;
    }

    if (x < (a + 1.0)){
        gser(&gamser, a, x, &gln, &doesNotConverge, &incorrectInput);

        if ( (doesNotConverge) || (incorrectInput) ){
            if (doesNotConverge)
                return DOES_NOT_CONVERGE;
            else 
                return INCORRECT_INPUT;
        }
        return gamser;
    } else {
        gcf(&gammcf, a, x, &gln, &doesNotConverge);

        if (doesNotConverge){
            return DOES_NOT_CONVERGE;
        }
        return 1.0 - gammcf;
    }
}


//gammq function
//returns negative values for errors
//basically the same as gammp, just 1-...
double gammq(double a, double x){
    double gamser; 
    double gammcf;
    double gln;
    bool doesNotConverge;
    bool incorrectInput;

    if ((x < 0.0) || (a <= 0.0)){
        //error - but this should never be reached because of wrapper
        return INCORRECT_INPUT;
    }

    if (x < (a + 1.0)){
        gser(&gamser, a, x, &gln, &doesNotConverge, &incorrectInput);

        if ( (doesNotConverge) || (incorrectInput) ){
            if (doesNotConverge)
                return DOES_NOT_CONVERGE;
            else 
                return INCORRECT_INPUT;
        }
        return 1.0 - gamser;
    } else {
        gcf(&gammcf, a, x, &gln, &doesNotConverge);
        if (doesNotConverge){
            return DOES_NOT_CONVERGE;
        }
        return gammcf;
    }
}
/*-------------------------------------------------------------------*/

//lower incomplete gamma function
//because gammp is scaled, this version "unscales it"
double gamma(double a, double x){
    double g = gammp(a, x);
    if (g > 0){
        return exp(gammaln(a)) * g;
    } else {
        return g;
    }
}

double gammaScaled(double x, double a){
    return gammp(a, x);
}

/*-------------------------------------------------------------------*/

//complementary error function
// 2/sqrt(pi) int_{x}^{infty} exp{-t^2} dt
//
// efficient implementation from Numerical Recipes in C
// Section 6.2
// Error is less than 1.2e-7
double erfc(double x){
    double t, z, ans;

    //using fabs instead of std::abs
    z = fabs(x);
    t = 1.0/(1.0 + 0.5*z);

    ans = t * exp(-z*z -1.26551223
            +t*(1.00002368
            +t*(0.37409196
            +t*(0.09678418
            +t*(-0.18628806
            +t*(0.27886807
            +t*(-1.13520398
            +t*(1.48851587
            +t*(-0.82215223
            +t*0.17087277)))))))));

    return x >= 0 ? ans : 2.0 - ans;
}

//later will check the speed
double erfc2(double x){
    double t, z, ans;

    //using fabs instead of std::abs
    z = fabs(x);
    t = 1.0/(1.0 + 0.5*z);

    const double cof[10] = {-1.26551223,
                           1.00002368,
                           0.37409196,
                           0.09678418,
                           -0.18628806,
                           0.27886807,
                           -1.13520398,
                           1.48851587,
                           -0.82215223,
                           0.17087277};

    ans = t * exp(-z*z + cof[0]  
            +t*(cof[1] 
            +t*(cof[2]
            +t*(cof[3]
            +t*(cof[4]
            +t*(cof[5]
            +t*(cof[6]
            +t*(cof[7]
            +t*(cof[8]
            +t*cof[9]
               )))))))));

    return x >= 0 ? ans : 2.0 - ans;
    //use a for loop rather
//     for (int i=9; i>0; i--){
//     
//     ans = t*exp(-z*z - cof[0])
//     }
}


//later will check the speed
double erfc3(double x){
    double t, z, ans;

    //using fabs instead of std::abs
    z = fabs(x);
    t = 1.0/(1.0 + 0.5*z);

    const double cof[10] = {-1.26551223,
                           1.00002368,
                           0.37409196,
                           0.09678418,
                           -0.18628806,
                           0.27886807,
                           -1.13520398,
                           1.48851587,
                           -0.82215223,
                           0.17087277};

    //use a for loop rather
    double aa = cof[9];
    for (int i=2; i<11; i--){
        aa *= t;
        aa += cof[10-i];
    }
    ans = t * exp(-z*z + aa);

    return x >= 0 ? ans : 2.0 - ans;
}


//error function
// 2/sqrt(pi) int_{0}^{x} exp{-t^2} dt
double erf(double x){
    return 1 - erfc(x);
}


//Normal cdf for any mu, sigma (check positive)
//NOTE: sigma, not sigma squared!
double normcdf(double x, double mu, double sigma){
    if (sigma <= 0){
        //error!
        return -1;
    }
    //double sqrt2 = 1.4142135623730951455;
    double sqrthalf = 0.70710678118654757274;

    double arg = sqrthalf * (x - mu) / sigma;
    return 0.5 * (1 + erf(arg));
}


//Beta distribution CONTINUED FRACTION
//Taken straight from Numerical Recipes in C, Section 6.4
//http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-4.pdf
double betacf(double a, double b, double x){
//     const int MAXIT=101;
//     const double EPS=3.0e-7;
//     const double FPMIN=1.0e-30;

    //initialising ints; m2 used for 2*m
    int m, m2;

    //initialising doubles
    double aa, c, d, del, h, qab, qam, qap;

    qab=a+b;
    qap=a+1.0;
    qam=a-1.0;
    c=1.0;
    d=1.0-qab*x/qap;

    if (fabs(d) < FPMIN) d=FPMIN;
    d=1.0/d;
    h=d;

    for (m=1; m < MAXIT; m++){
        m2=2*m;

        aa=m*(b-m)*x/((qam+m2) * (a+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c = 1.0+aa/c; 
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        h *= d*c;

        aa = -(a+m)*(qab+m)*x/((a+m2) * (qap+m2));
        d=1.0+aa*d;
        if (fabs(d) < FPMIN) d=FPMIN;
        c = 1.0+aa/c; 
        if (fabs(c) < FPMIN) c=FPMIN;
        d=1.0/d;
        del=d*c;
        h *= del;

        if (fabs(del-1.0) < EPS) break;
    }
    if (m > MAXIT)
        return(DOES_NOT_CONVERGE);
    return h;
}


//Incomplete beta function
//i.e. computes 
//$$
//I_{x}(a,b) = \frac{1}{B(a, b)} \int_{0}^{x} t^{a-1} (1-t)^{b-1} dt
//(a, b > 0)
//$$
//
//Taken straight from Numerical Recipes in C, Section 6.4
//http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-4.pdf
double betai(double a, double b, double x){
    double bt;
    double cf;

    if (x < 0.0 || x > 1.0 || a <= 0 || b <= 0)
        return(INCORRECT_INPUT);
    if (x == 0.0 || x == 1.0) 
        bt=0.0;
    else
        bt=exp(gammaln(a+b)-gammaln(a)-gammaln(b) + a*log(x) + b*log(1.0-x));

    if(x < (a+1.0)/(a+b+2.0)){
        cf = betacf(a,b,x);
        if (cf==DOES_NOT_CONVERGE){
            return DOES_NOT_CONVERGE;
        }
        return bt*cf/a;
    
    } else{
        cf = betacf(b, a, 1.0-x);
        if (cf==DOES_NOT_CONVERGE){
            return DOES_NOT_CONVERGE;
        }
        return 1.0-bt*cf/b;
    }
}

//Computes the cdf of Student's t-distribution by
//using the incomplete beta function
//See Numerical Recipes in C, Section 6.4
//http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-4.pdf
//*NB*: Found an error in the above link; in Equation (6.4.9)
//      there should be a factor of 0.5 in front of the I function
double tcdf(double t, double nu){
    //nu needs to be strictly larger than 0
    if (nu <= 0){
        return INCORRECT_INPUT;
    }

    //following relation between t-dist cdf and betai
    double val = betai( nu*0.5, 0.5, nu/(nu+t*t) );

    //need to check for bad output, incorrect input/does not converge
    if (val < 0)
        return val;

    //otherwise, return correct answer
    //*NB*: The extra factor of 0.5 is necessary! Typo in Num. Recipes in C!
    return 1.0 -  0.5*val;
}

#endif
