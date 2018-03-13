
##' auxiliary function for applying the function 'IntegrandInterpolationMethod', approximate solution of the filter integral by linear interpolation
##' 
##' @title IntegrateLinearInterpolation
##' @param tx.pos time points of the time series covered by the filter window and used during integrand interpolation
##' @param x.pos product of filter weight and observation at time points of the time series covered by the filter window and used during integrand interpolation
##' @return estimate of the filter integral
##' @examples 
##' 
##' @author
##' @export


IntegrateLinearInterpolation <- function(tx.pos, x.pos)
{
	diff.tx <- diff(tx.pos)
	diff.x <- abs(diff(x.pos))
	A <- numeric()
	#
	for(n in 1:length(diff.tx)){
		#
		tmpA <- diff.tx[n] * ((x.pos[n] + x.pos[n+1])/2)
		#
		A[n] <- tmpA
	}
	return(sum(A))
}
