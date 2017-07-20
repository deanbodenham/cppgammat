context("Test5: Student's t cdf")

tol <<- 1.2e-7
#tol <<- 1.2e-8

r_tcdf <- function(t, nu){ return(pt(q=t, df=nu)) }

#r_ibeta <- function(a,b,x){ pbeta(q=x,shape1=a,shape2=b) }
r_tcdf2 <- function(t, nu){ 
    x <- nu/(nu+t^2)
    a <- nu/2
    b <- 0.5
    val <- 1.0 - 0.5 * pbeta(q=x, shape1=a, shape2=b)
    return(val)
}
    


test_that("Testing alternative formula for Student's t cdf", {
    tVec <- c(1.2, 3.12, 6.43, 8.752)
    nuVec <- c(0.5, 2.32, 4.5, 6.975)

    for (i in seq_along(tVec)){
        t <- tVec[i]
        nu <- nuVec[i]
        expect_equal(r_tcdf2(t, nu), r_tcdf(t, nu), tol=tol)
#        print(r_tcdf(t, nu))
#        print(r_tcdf2(t, nu))
#        print(tcdfRcpp(t, nu))
    }
})


test_that("Testing Student's t cdf", {
    tVec <- c(1.2)
    nuVec <- c(3.0)

    for (i in seq_along(tVec)){
        t <- tVec[i]
        nu <- nuVec[i]
        expect_equal(tcdfRcpp(t, nu), r_tcdf(t, nu), tol=tol)
    }

})


test_that("Testing for bad values 1", {
    t <- 0.3
    nu <- -1
    expect_error(tcdfRcpp(t, nu), "strictly greater than 0")
})


test_that("Testing for bad values 1", {
    t <- 0.3
    nu <- 0
    expect_error(tcdfRcpp(t, nu), "strictly greater than 0")
})
