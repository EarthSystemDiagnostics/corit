##' auxiliary function to get the filter weights of the Gaussian filter of a certain cut-off frequency
##' 
##' @title gauss.weights
##' @param fc cut-off frequency
##' @return filter weights 
##' @examples 
##' 
##' @author
##' @export


gauss.weights <- function(fc)
{
	h <- 1/ (fc * pi * sqrt(2/log(2)))
	filt.position <- seq(from=round(-3*h), to=round(3*h), by=1)
	filter <- 1/sqrt(2 * pi * h^2) * (exp(-((filt.position^2))/(2 * h^2)))
	#
	return(filter/sum(filter))
}
