##' auxiliary function for applying the functions 'GenNullHypPair' & 'CorQuantilesNullHyp', detrending a time series by removing a linear trend
##' 
##' @title detrTimser
##' @param timeseries time series to be detrendend
##' @param detr TRUE (detrending) or FALSE (no detrending)
##' @return zoo-object of the (not) detrendend time series 
##' @examples 
##' 
##' @author
##' @export


detrTimser <- function(timeseries, detr=TRUE)
{
	if(detr==TRUE){
		result <- zoo(lm(coredata(timeseries)~index(timeseries), na.action=na.omit)$residuals, order.by=index(timeseries))
	}else{
		result <- timeseries
	}
	#
	return(result)
}
