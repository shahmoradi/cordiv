# This is a variant of the function corrplot.mixed in R package 'corrplot'. It makes some changes to the color and format of the correlation matrix figure, that I need for the manuscript.
# Amir Shahmoradi, Wednesday 6:54 PM, Nov 19 2014, iCMB, UT Austin
#
mycpm <- function (corr, lower = "number", upper = "circle", tl.pos = c("d", 
    "lt", "n"), diag = c("n", "l", "u"), bg = "white", addgrid.col = "gray", 
    ...) 
{
    diag <- match.arg(diag)
    tl.pos <- match.arg(tl.pos)
    n <- nrow(corr)
    corrplot(corr, type = "upper", method = lower, diag = TRUE, 
        tl.pos = tl.pos, ...)
    corrplot(corr, add = TRUE, type = "lower", method = upper, 
        diag = (diag == "l"), tl.pos = "n", cl.pos = "n", col = 'black', ...)
    if (diag == "n" & tl.pos != "d") {
        symbols(1:n, n:1, add = TRUE, bg = bg, fg = addgrid.col, 
            inches = FALSE, squares = rep(1, n))
    }
}