context("Test 3: Normal cdf and error function")

tol <<- 1.2e-7

erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

test_that("Testing erf for positive values", {
    xVec <- c(1.2, 3.4, 5.6)

    for (x in xVec){
        expect_equal(erfRcpp(x), erf(x), tol=tol)
    }
})


test_that("Testing erf for negative values", {
    xVec <- c(-1.2, -3.4, -5.6)

    for (x in xVec){
        expect_equal(erfRcpp(x), erf(x), tol=tol)
    }
})


#actually, tol is only for the error function;
#  hopefully it carries over to the normal cdf
test_that("Testing normal cdf", {
    xVec <- c(1.3, -7.3, 4.3, 9.222)
    muVec <- c(0, 2.3, 5.4, -2)
    sigmaVec <- c(1, 1, 2, 3)

    for (i in seq_along(xVec)){
        x <- xVec[i]
        mu <- muVec[i]
        sigma <- sigmaVec[i]

        expect_equal(normcdfRcpp(x, mu, sigma), pnorm(x, mu, sigma), tol=tol)
    }

})


