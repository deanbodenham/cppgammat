# C++ implementation of Gamma and Beta Functions


Testing Equation (6.1.36) in Abramowitz and Stegun, page 257 for Gamma function:
[http://people.math.sfu.ca/~cbm/aands/page_257.htm](http://people.math.sfu.ca/~cbm/aands/page_257.htm)

Which is used after scaling `x` in the CDF of the Chi-squared distribution in Equation (26.4.1):
[http://people.math.sfu.ca/~cbm/aands/page_940.htm](http://people.math.sfu.ca/~cbm/aands/page_940.htm)

Also implemented incomplete Beta function and Student's *t*-distribution CDF.

These implementations were obtained from:  
[Numerical Recipes in C Section 6.1](http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-1.pdf)  
[Numerical Recipes in C Section 6.2](http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-2.pdf)  
[Numerical Recipes in C Section 6.4](http://www.aip.de/groups/soe/local/numres/bookcpdf/c6-4.pdf)

**Note:** There is a small error in Equation (6.4.9), there should be an extra factor of `0.5` in front of the incomplete Beta function, i.e.
```
1.0 - 0.5 * I(nu/(nu+t^2), nu/2, 1/2)

```



## Tests

Part of the goal of this project was for it to be an exercise to run tests in both R and C++, using the `testthat` and `googletest` frameworks, respectively.
