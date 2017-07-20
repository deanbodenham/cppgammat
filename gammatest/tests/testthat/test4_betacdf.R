context("Test 4: Beta cdf")

tol <<- 1.2e-9

#r_beta_cdf <- function(a, b, x){
#    return( pbeta(q=x, shape1=b, shape2=2*a) )
#    return( pbeta(q=x, shape1=a, shape2=b) )
#}

r_ibeta <- function(a,b,x){ pbeta(q=x,shape1=a,shape2=b) }
#r_ibeta <- function(a,b,x){ pbeta(q=x,shape1=a,shape2=b)*beta(a,b) }
#r_ibeta <- function(a,b,x){ pbeta(q=x,shape1=a,shape2=b)*beta(a,b) }
#r_ibeta <- function(a,b,x){ pbeta(x,a,b)*beta(a,b) }

test_that("Testing Beta cdf", {
    aVec <- c(1, 2, 3, 3) 
    bVec <- c(2, 1, 4, 4)
    xVec <- c(0.5, 0.6, 0.7, 0.5)

    for (i in seq_along(xVec)){
        a <- aVec[i]
        b <- bVec[i]
        x <- xVec[i]
        expect_equal( betaiRcpp(a, b, x), r_ibeta(a, b, x), tol=tol )
    }
})



test_that("Testing Beta cdf 2", {
    aVec <- c(3.12) 
    bVec <- c(4.765)
    xVec <- c(0.5678)

    for (i in seq_along(xVec)){
        a <- aVec[i]
        b <- bVec[i]
        x <- xVec[i]
        expect_equal( betaiRcpp(a, b, x), r_ibeta(a, b, x), tol=tol )
    }
})


test_that("Testing Beta cdf limits", {
    aVec <- c(3.12, 1) 
    bVec <- c(4.765, 2)
    xVec <- c(0, 1)

    for (i in seq_along(xVec)){
        a <- aVec[i]
        b <- bVec[i]
        x <- xVec[i]
        expect_equal( betaiRcpp(a, b, x), r_ibeta(a, b, x), tol=tol)
    }
#    print(r_ibeta(4, 5, 1))
#    print(r_ibeta(4, 5, 0))
})


test_that("Testing Beta function for bad x values", {
    aVec <- c(2, 2, -0.9, 1) 
    bVec <- c(3, 3, 1, -0.9)
    xVec <- c(1.1, -0.9, 0.5, 0.5) 

    for (i in seq_along(xVec)){
        a <- aVec[i]
        b <- bVec[i]
        x <- xVec[i]
        expect_error(betaiRcpp(a, b, x), "both be strictly")
    }
})
