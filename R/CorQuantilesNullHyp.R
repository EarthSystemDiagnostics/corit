
##' Monte Carlo based procedure to estimate the quantiles of the correlation of independent noise which belong to the methods and applied filter to 
##' estimate the time scale dependent correlation of irregularly sampled time series
##' 
##' @title CorQuantilesNullHyp
##' @param timser1,timser2 time series (zoo-objects) for which the Null Hypothesis of the correlation is to be applied
##' @param beta.noise1,beta.noise2 fitted slopes from the power spectrum of timser1 and timser2
##' @param detr TRUE for removing a linear trend, else FALSE
##' @param rep number of repetitions during the Monte Carlo procedure
##' @param quant quantiles to be estimated
##' @param method method to handle irregularity of sampling (DirectFiltering, IntegrandInterpolationMethod, InterpolationMethod)
##' @param appliedFilter time domain filter (gauss, runmean, lowpass)
##' @param fc cut-off frequency of the applied filter
##' @param tn output vector (time) of the filtered data (only used in case of DirectFiltering and IntegrandInterpolationMethod)
##' @param dt regular inter-observation time step of the interpolation (only used in case of InterpolationMethod)
##' @param int.method kind of interpolation (linear, nearest neighbor) (only used in case of InterpolationMethod)
##' @param k scaling factor to define the sharpness of the lowpass
##' @return $corPair: estimated correlations during the Monte Carlo procedure, $Quantile: quantiles to be estimated
##' @examples 
##' timeseries1 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timeseries2 <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' slopes <- estimateTimserSlopes(timeseries1,timeseries2,1)
##' CorQuantilesNullHyp(timser1=timeseries1, timser2=timeseries2, beta.noise1=slopes$s1, beta.noise2=slopes$s2, detr=FALSE, rep=1000, quant=c(0.05,0.95), method="InterpolationMethod", appliedFilter="gauss", fc=1/200, tn=NA, dt=10, int.method="linear", k=NA)
##' @author
##' @export



CorQuantilesNullHyp <- function(timser1, timser2, beta.noise1, beta.noise2, detr, rep, quant, method=c("InterpolationMethod","DirectFiltering","IntegrandInterpolationMethod"), appliedFilter=c("gauss", "runmean", "lowpass"), fc, tn=seq(from=10,to=max(c(index(timser1),index(timser2))),by=10), dt, int.method=c("linear","nearest"), k=5)
{
	n <- max(c(max(index(timser1)), max(index(timser2))))
	time1 <- index(timser1)
	time2 <- index(timser2)
	#
	corNullHypPair <- QuantCorPair <- list()
	tmpCorNullHypPair <- numeric(length=rep)
	
	for(p in 1:length(fc)){
		#print(paste("p=",p,sep=""))
		for(i in 1:rep){
			#print(i)
			tmp <- GenNullHypPair(n, beta.noise1, beta.noise2, time1, time2, detr)
			#
			if(method=="InterpolationMethod"){
				tmpCorNullHypPair[i] <- cor.test(InterpolationMethod(tmp$y1,fc[p],dt,n,int.method,appliedFilter,k), InterpolationMethod(tmp$y2,fc[p],dt,n,int.method,appliedFilter,k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
			}
			if(method=="DirectFiltering"){
				tmpCorNullHypPair[i] <- cor.test(DirectFiltering(tmp$y1,fc[p],tn,appliedFilter,k), DirectFiltering(tmp$y2,fc[p],tn,appliedFilter,k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
			}
			if(method=="IntegrandInterpolationMethod"){
				tmpCorNullHypPair[i] <- cor.test(IntegrandInterpolationMethod(tmp$y1,fc[p],tn,appliedFilter,k), IntegrandInterpolationMethod(tmp$y2,fc[p],tn,appliedFilter,k), method="pearson", alternative="two.sided", na.action=TRUE)$estimate
			}
		}
		corNullHypPair[[p]] <- tmpCorNullHypPair
		QuantCorPair[[p]] <- quantile(corNullHypPair[[p]],quant)
	}
	#
	return(list(corPair=corNullHypPair, Quantile=QuantCorPair))
}
