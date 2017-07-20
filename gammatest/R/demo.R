#' Demo of basic timesTwo function
#' 
#' A simple driver for the timesTwoRcpp function. No parameters.
#'
#' @examples
#' demoTimesTwoRcpp()
#'
#' @export
demoTimesTwoRcpp <- function(){
    xVec <- c(4, 2.3)
    for (i in seq_along(xVec)){
        x <- xVec[i]
        x2 <- timesTwoRcpp(x);
        cat("e.g. ", i, ": ", x, " times two is: ", x2, "\n", sep="")
    }
}
