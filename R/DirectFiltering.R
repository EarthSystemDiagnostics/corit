
##' direct filtering with application of the time domain filter to irregular time series at regular filter positions
##' 
##' @title DirectFiltering
##' @param X time series (zoo-object)
##' @param fc cut-off frequency
##' @param tn output vector (time) of the filtered data
##' @param appliedFilter time domain filter (gauss, runmean, lowpass)
##' @param k scaling factor to define the sharpness of the lowpass
##' @return filtered time series (zoo-object)
##' @examples 
##' timser.irreg <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timser.irregSmooth <- DirectFiltering(X=timser.irreg, fc=1/50, tn=seq(from=10,to=1000,by=10), appliedFilter="lowpass", k=5)
##' plot(timser.irreg,type="l",col="black")
##' lines(timser.irregSmooth,col="limegreen",lwd=3)
##' @author
##' @export



DirectFiltering <- function(X, fc, tn, appliedFilter=c("gauss","runmean","lowpass"), k=5)
{
	if (!is.zoo(X)){stop("time series must be zoo-object")}
	#
	tx <- index(X)
	x <- coredata(X)
	Xsmooth <- rep(NA,length(tn))
	#
		if(appliedFilter=="gauss"){
			h <- 1/ (fc*pi*sqrt(2/log(2)))
			#
			for(i in 1:length(tn)){	
			if((tn[i] >= (min(tx) + 3*h)) & (tn[i] <= (max(tx) - 3*h)))
				{
					filter <- 1/sqrt(2 * pi * h^2) * (exp(-((tx - tn[i])^2)/(2 * h^2)))
					Xsmooth[i] <- sum(filter * x)/ sum(filter)
				}
			}
		}
		if(appliedFilter=="runmean"){
			filt.length <- 1/(2*fc)
			if (filt.length <= max(diff(tx))){warning("filt.length must be larger than max time step. output contains NAs.")}
			#
			for(i in 1:length(tn)){
				if((tn[i] > (min(tx) + (filt.length/2))) & (tn[i] < (max(tx) - (filt.length/2))))
				{
					ind.elems <- which((tx >= (tn[i] - (filt.length/2))) & (tx <= (tn[i] + (filt.length/2))))
					Xsmooth[i] <- mean(x[ind.elems], na.rm=TRUE)
				}
			}
		}
		if(appliedFilter=="lowpass"){	
			filt.length <- k/fc
				if((filt.length %% 2)==0){filt.length <- filt.length+1}
				if(filt.length <= max(diff(tx))){warning("filt.length must be larger than max time step. output contains NAs.")}
			if((filt.length %% 2)==0){filt.length <- filt.length+1}
			LP.w <- Lowpass(omega.c=fc, n=filt.length, sample=1)
			LP.diff <- seq(from= -(filt.length-1)/2, to=(filt.length-1)/2, by=1)
			#
			for(i in 1:length(tn)){
				if((tn[i] >= (min(tx) + ((filt.length-1)/2))) & (tn[i] <= (max(tx) - ((filt.length-1)/2))))
				{
					diff.i<-tx-tn[i]
					filt.ind <- which((diff.i >= min(LP.diff)) & (diff.i <= max(LP.diff)))
					if(length(filt.ind)==0)
					{
						Xsmooth[i] <- NA
					} 
					else 
					{
						filter <- approx(LP.diff, LP.w, xout=diff.i[filt.ind])$y
						Xsmooth[i] <- sum(filter * x[filt.ind])/ sum(filter)
					}
				}
			}
		}
	#
	return(zoo(Xsmooth, order.by=tn))
}



