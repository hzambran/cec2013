cec2013 <- function (i, x) {
    if (is.numeric(i) && i >= 1 && i <= 28) {
        if (is.vector(x)) {
            row <- 1; col <- length(x)
        } else if (is.matrix(x)) {
            row <- nrow(x); col <- ncol(x)
        } else {
            stop("x should be a vector or a matrix")
        }
        if (!(col %in% c(2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100))) {
            stop("only 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90 and 100 variables are supported")
        }
        extdatadir <- system.file("extdata", package = "cec2013")
        f <- .C("cec2013", extdatadir = as.character(extdatadir), 
                i = as.integer(i), x = as.double(x), row = as.integer(row),
                col = as.integer(col), f = double(row), 
                PACKAGE = "cec2013")$f
    } else {
        stop("i should be an integer between 1 and 28")
    }
    return(f)
}
