
##' integrand interpolation of the filter integral during filtering irregular time series at regular positions
##' 
##' @title IntegrandInterpolationMethod
##' @param X time series (zoo-object)
##' @param fc cut-off frequency
##' @param tn output vector (time) of the filtered data
##' @param appliedFilter time domain filter (gauss, runmean, lowpass)
##' @param k scaling factor to define the sharpness of the lowpass
##' @return filtered time series (zoo-object)
##' @examples 
##' timser.irreg <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timser.irregSmooth <- IntegrandInterpolationMethod(X=timser.irreg, fc=1/50, tn=seq(from=10,to=1000,by=10), appliedFilter="lowpass", k=5)
##' plot(timser.irreg,type="l",col="black")
##' lines(timser.irregSmooth,col="limegreen",lwd=3)
##' @author
##' @export



IntegrandInterpolationMethod <- function(X, fc, tn, appliedFilter=c("gauss","runmean","lowpass"), k=5)
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
					filt.pos <- which((tx >= (tn[i]-3*h)) & (tx <= (tn[i]+3*h)))
					filter <- 1/sqrt(2 * pi * h^2) * (exp(-((tx[filt.pos] - tn[i])^2)/(2 * h^2)))
					filter.product <- (filter * x[filt.pos])
					#
					if((tn[i] - 3*h) == min(tx))
						{
							tx.pos <- c(tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(filter.product, 0)
						}
					if((tn[i] + 3*h) == max(tx))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos])
							x.pos <- c(0, filter.product)
						}
					if(((tn[i] - 3*h) != min(tx)) & ((tn[i] + 3*h) != max(tx)))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(0, filter.product, 0)
						}
					#
					Xsmooth[i] <- IntegrateLinearInterpolation(tx.pos, x.pos)
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
					filt.pos <- which((tx >= (tn[i]-(filt.length/2))) & (tx <= (tn[i]+(filt.length/2))))
					filter.product <- x[filt.pos] * (1/filt.length)
					#
					if((tn[i] - filt.length/2) == min(tx))
						{
							tx.pos <- c(tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(filter.product, 0)
						}
					if((tn[i] + filt.length/2) == max(tx))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos])
							x.pos <- c(0, filter.product)
						}
					if(((tn[i] - filt.length/2) != min(tx)) & ((tn[i] + filt.length/2) != max(tx)))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(0, filter.product, 0)
						}
					#
					Xsmooth[i] <- IntegrateLinearInterpolation(tx.pos, x.pos)
				}
			}	
		}
		if(appliedFilter=="lowpass"){
			filt.length <- k/fc
				if((filt.length %% 2)==0){filt.length <- filt.length+1}
				if(filt.length <= max(diff(tx))){warning("filt.length must be larger than max time step. output contains NAs.")}
			if((filt.length %% 2)==0){filt.length <- filt.length+1}
			delta <- 4*pi/ filt.length
			#
			for(i in 1:length(tn)){
				if((tn[i] >= (min(tx) + ((filt.length-1)/2))) & (tn[i] <= (max(tx) - ((filt.length-1)/2))))
				{
					filt.pos <- which((tx >= (tn[i]-(filt.length-1)/2)) & (tx <= (tn[i]+(filt.length-1)/2)))
					filter <- (sin((tx[filt.pos]-tn[i])*fc*2*pi)/ ((tx[filt.pos]-tn[i])*pi)) * (sin((tx[filt.pos]-tn[i])*delta/2)/ ((tx[filt.pos]-tn[i])*delta/2))
					filter.product <- (filter * x[filt.pos])
					#
					if((tn[i] - (filt.length-1)/2) == min(tx))
						{
							tx.pos <- c(tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(filter.product, 0)
						}
					if((tn[i] + (filt.length-1)/2) == max(tx))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos])
							x.pos <- c(0, filter.product)
						}
					if(((tn[i] - (filt.length-1)/2) != min(tx)) & ((tn[i] + (filt.length-1)/2) != max(tx)))
						{
							tx.pos <- c(tx[min(filt.pos)-1], tx[filt.pos], tx[max(filt.pos)+1])
							x.pos <- c(0, filter.product, 0)
						}
					#
					Xsmooth[i] <- IntegrateLinearInterpolation(tx.pos, x.pos)
				}
			}
		}
	#
	return(zoo(Xsmooth, order.by=tn))	
}

