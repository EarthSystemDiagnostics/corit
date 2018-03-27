##' generating a pair of subsampled power-law time series 
##' 
##' @title GenNullHypPair
##' @param n length of the power-law time series which have to be generated
##' @param beta.noise1,beta.noise2 slopes from the power spectrum of the time series to be generated
##' @param time1,time2 time points at which the simulated power-law has to be subsampled
##' @param detr TRUE for removing a linear trend, else FALSE
##' @return zoo-objects of the time series 
##' @examples 
##' time.series <- GenNullHypPair(1000, 0.5, 0.4, sort(runif(100, min=1, max=1000)), sort(runif(100, min=1, max=1000)), FALSE)
##' @author
##' @export


GenNullHypPair <- function(n, beta.noise1, beta.noise2, time1, time2, detr)
{
	regtimeseries1 <- SimPowerlaw(beta.noise1,n)
	regtimeseries2 <- SimPowerlaw(beta.noise2,n)
	#
	timeseries1 <- SubsampleTimeseriesBlock(ts=regtimeseries1, timepoints=time1)
	timeseries2 <- SubsampleTimeseriesBlock(ts=regtimeseries2, timepoints=time2)
	#
	return(list(y1=detrTimser(zoo(timeseries1,order.by=time1), detr), y2=detrTimser(zoo(timeseries2,order.by=time2), detr)))
}
