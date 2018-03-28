
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
##' @param filt.output TRUE for returning correlation ($cor) and filter results ($ft1, $ft2), FALSE for only returning correlation
##' @examples 
##' timeseries1 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timeseries2 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' CorIrregTimser(timser1=timeseries1, timser2=timeseries2, detr=FALSE, method="InterpolationMethod", appliedFilter="gauss", fc=1/200, tn=NA, dt=10, int.method="linear", k=NA, filt.output=TRUE)
##' @author
##' @export



CorIrregTimser <- function(timser1, timser2, detr, method=c("InterpolationMethod","DirectFiltering","IntegrandInterpolationMethod"), appliedFilter=c("gauss", "runmean", "lowpass"), fc, tn=seq(from=10,to=max(c(index(timser1),index(timser2))),by=10), dt, int.method=c("linear","nearest"), k=5, filt.output)
{
	n <- max(c(max(index(timser1)), max(index(timser2))))
	#
	if(method=="InterpolationMethod"){
		filtTimser1 <- InterpolationMethod(detrTimser(timser1, detr), fc, dt, n, int.method, appliedFilter, k)
		filtTimser2 <- InterpolationMethod(detrTimser(timser2, detr), fc, dt, n, int.method, appliedFilter, k)
		res <- cor.test(filtTimser1, filtTimser2, method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	if(method=="DirectFiltering"){
		filtTimser1 <- DirectFiltering(detrTimser(timser1, detr), fc, tn, appliedFilter, k)
		filtTimser2 <- DirectFiltering(detrTimser(timser2, detr), fc, tn, appliedFilter, k)
		res <- cor.test(filtTimser1, filtTimser2, method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	if(method=="IntegrandInterpolationMethod"){
		filtTimser1 <- IntegrandInterpolationMethod(detrTimser(timser1, detr), fc, tn, appliedFilter, k)
		filtTimser2 <- IntegrandInterpolationMethod(detrTimser(timser2, detr), fc, tn, appliedFilter, k)
		res <- cor.test(filtTimser1, filtTimser2, method="pearson", alternative="two.sided", na.action=TRUE)$estimate
	}
	#
	if(filt.output==TRUE){return(list(cor = res, ft1 = filtTimser1, ft2 = filtTimser2))}
	if(filt.output==FALSE){return(res)}
}
