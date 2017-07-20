context("Test 1: Test timesTwoRcpp")

test_that("integer timesTwoRcpp", {
    expect_equal(timesTwoRcpp(1), 2)
    expect_equal(timesTwoRcpp(1234), 2468)
})


test_that("double timesTwoRcpp", {
    expect_equal(timesTwoRcpp(3.421), 6.842)
})


test_that("negative timesTwoRcpp", {
    expect_equal(timesTwoRcpp(-1.3056), -2.6112)
})


test_that("up to precision", {
    expect_equal(timesTwoRcpp(1.234), 2.468001, tol=1e-6)
})
