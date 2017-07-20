context("Test 2: The gamma function")

incompleteGamma <- function(a, x){
    return(pgamma(x, a) * gamma(a))
}


test_that("gammaln, various", {
    xVec <- c(1, 2, 2, 3, 5, 4, 3.2, 1.3, 10.72)
    for (x in xVec){
        expect_equal(gammalnRcpp(x), log(gamma(x)))
    }
})


test_that("Incomplete gamma function Scaled", {
    tol <- 3.0e-7
    xVec <- c(1, 2, 3, 2, 4.1, 5, 5)
    aVec <- c(1, 2, 2, 3, 5, 4.1, 3.9)
    for (i in seq_along(xVec)){
        x <- xVec[i]
        a <- aVec[i]
#        cat("\n", "a: ", a, ", x: ", x, "\n\n")
        expect_equal(gammaScaledRcpp(x, a), pgamma(x, a), tol=tol)
    }

})

test_that("Incomplete gamma function", {
    tol <- 3.0e-7
    xVec <- c(1, 2, 3, 2, 4, 5)
    aVec <- c(1, 2, 2, 3, 5, 4)
    for (i in seq_along(xVec)){
        x <- xVec[i]
        a <- aVec[i]
#        cat("\n", "a: ", a, ", x: ", x, "\n\n")
        expect_equal(gammaRcpp(a, x), incompleteGamma(a, x), tol=tol)
    }

})

