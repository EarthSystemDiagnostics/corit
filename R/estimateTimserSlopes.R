
##' estimating the slope from a power spectrum of a pair of time series by a linear fit in the frequency range of 
##' the inverse half of the length of the overlapping time window of both time series and the inverse twofold maximum
##' mean resolution of both time series
##' 
##' @title estimateTimserSlopes
##' @param timeseries1,timeseries2 (zoo-object)
##' @param int.step equidistant time step of the interpolated irregular time series
##' @return $s1,$s2 fitted slopes of both time series 
##' @examples 
##' timeseries1 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timeseries2 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' estimateTimserSlopes(timeseries1,timeseries2,1)
##' @author
##' @export



estimateTimserSlopes <- function(timeseries1, timeseries2, int.step)
{
	tmp1 <- approx(x=index(timeseries1), y=coredata(timeseries1), xout=seq(from=round(min(index(timeseries1)),0), to=round(max(index(timeseries1)),0), by=int.step), method="linear")
	approxTimser1 <- zoo(tmp1$y, order.by=tmp1$x)
	specTimser1 <- spectrum(approxTimser1[which(!is.na(approxTimser1))],plot=FALSE)
	#
	tmp2 <- approx(x=index(timeseries2), y=coredata(timeseries2), xout=seq(from=round(min(index(timeseries2)),0), to=round(max(index(timeseries2)),0), by=int.step), method="linear")
	approxTimser2 <- zoo(tmp2$y, order.by=tmp2$x)
	specTimser2 <- spectrum(approxTimser2[which(!is.na(approxTimser2))],plot=FALSE)
	#
	slope1 <- SlopeFit(specTimser1, freq.start=1/(0.5*diff(c(max(min(index(approxTimser1)),min(index(approxTimser2))),min(max(index(approxTimser1)),max(index(approxTimser2)))))), freq.end=1/max(c((2*mean(diff(index(timeseries1)))),(2*mean(diff(index(timeseries2)))))),bDebug=FALSE,i.fStart=1)
	#
	slope2 <- SlopeFit(specTimser2, freq.start=1/(0.5*diff(c(max(min(index(approxTimser1)),min(index(approxTimser2))),min(max(index(approxTimser1)),max(index(approxTimser2)))))), freq.end=1/max(c((2*mean(diff(index(timeseries1)))),(2*mean(diff(index(timeseries2)))))),bDebug=FALSE,i.fStart=1)
	#
	return(list(s1=abs(slope1$slope), s2=abs(slope2$slope)))
}
