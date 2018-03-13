##' auxiliary function for applying the function 'InterpolationMethod', approximate a time series using the nearest neighbour
##' 
##' @title approx.nearest
##' @param x numeric vector giving the coordinates of the points to be interpolated
##' @param y corresponding values to x
##' @param xi set of numeric values specifying where interpolation has to take place
##' @return zoo-object containing the interpolated data
##' @examples 
##' 
##' @author
##' @export


approx.nearest <- function(x, y, xi)
{
	result <- list()
	result$x <- xi
	result$y <- approx(c(x[1], x+c(diff(x)/2, 0)), c(y[1], y), xi, method="constant", f=1)$y
	output <- zoo(result$y, order.by=result$x)
	#
	return(output)
}
