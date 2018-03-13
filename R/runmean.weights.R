##' auxiliary function to get the filter weights of the moving average filter of a certain cut-off frequency
##' 
##' @title runmean.weights
##' @param fc cut-off frequency
##' @return filter weights
##' @examples 
##' 
##' @author
##' @export


runmean.weights <- function(fc)
{
	filt.length <- 1/(2*fc)
		if((filt.length %% 2)==0){filt.length <- filt.length+1}
	filter <- rep.int(1, filt.length)
	#
	return(filter/sum(filter))
}
