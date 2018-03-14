
##' estimating the time scale dependent correlation of irregularly sampled time series
##' 
##' @title CorIrregTimser
##' @param timser1,timser2 time series (zoo-objects) 
##' @param detr TRUE for removing a linear trend, else FALSE
##' @param method method to handle irregularity of sampling (DirectFiltering, IntegrandInterpolationMethod, InterpolationMethod)
##' @param appliedFilter time domain filter (gauss, runmean, lowpass)
##' @param fc cut-off frequency of the applied filter
##' @param tn output vector (time) of the filtered data (only used in case of DirectFiltering and IntegrandInterpolationMethod)
##' @param dt regular inter-observation time step of the interpolation (only used in case of InterpolationMethod)
##' @param int.method kind of interpolation (linear, nearest neighbor) (only used in case of InterpolationMethod)
##' @param k scaling factor to define the sharpness of the lowpass
##' @examples 
##' timeseries1 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timeseries2 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' CorIrregTimser(timeseries1, timeseries2, FALSE, "InterpolationMethod", "gauss", 1/200, NA, 10, "linear", NA)
##' @author
##' @export



CorIrregTimser <- function(timser1, timser2, detr, method=c("InterpolationMethod","DirectFiltering","IntegrandInterpolationMethod"), appliedFilter=c("gauss", "runmean", "lowpass"), fc, tn=seq(from=10,to=max(c(index(timser1),index(timser2))),by=10), dt, int.method=c("linear","nearest"), k=5)
{
	n <- max(c(max(index(timser1)), max(index(timser2))))
	#
	if(method=="InterpolationMethod"){
		res <- cor.test(InterpolationMethod(detrTimser(timser1, detr), fc, dt, n, int.method, appliedFilter, k), InterpolationMethod(detrTimser(timser2, detr), fc, dt, n, int.method, appliedFilter, k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	if(method=="DirectFiltering"){
		res <- cor.test(DirectFiltering(detrTimser(timser1, detr), fc, tn, appliedFilter, k), DirectFiltering(detrTimser(timser2, detr), fc, tn, appliedFilter, k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	if(method=="IntegrandInterpolationMethod"){
		res <- cor.test(IntegrandInterpolationMethod(detrTimser(timser1, detr), fc, tn, appliedFilter, k), IntegrandInterpolationMethod(detrTimser(timser2, detr), fc, tn, appliedFilter, k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	#
	return(res)
}
