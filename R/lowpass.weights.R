##' auxiliary function to get the filter weights of the low-pass filter (based on Bloomfield, 1976) of a certain cut-off frequency, necessitates the auxiliary function 'Lowpass'
##' 
##' @title lowpass.weights
##' @param fc cut-off frequency
##' @param k scaling factor to define the sharpness of the lowpass
##' @return filter weights
##' @examples 
##' 
##' @author
##' @export


lowpass.weights <- function(fc, k=5)
{
	filt.length <- k/fc
		if((filt.length %% 2)==0){filt.length <- filt.length+1}
	LP.w <- Lowpass(omega.c=fc, n=filt.length, sample=1)
	#
	return(LP.w/sum(LP.w))
}
