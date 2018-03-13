
##' linear or nearest neighbor interpolation of the irregular time series before convolution with time domain filter
##' 
##' @title InterpolationMethod
##' @param X time series (zoo-object)
##' @param fc cut-off frequency
##' @param dt regular inter-observation time step of the interpolation
##' @param timser.length oldest age of the time series
##' @param int.method kind of interpolation (linear, nearest neighbor)
##' @param appliedFilter time domain filter (gauss, runmean, lowpass)
##' @param k scaling factor to define the sharpness of the lowpass
##' @return filtered time series (zoo-object)
##' @examples 
##' timser.irreg <- zoo(rnorm(100), order.by=sort(runif(100,min=1,max=1000)))
##' timser.irregSmooth <- InterpolationMethod(X=timser.irreg, fc=1/50, dt=10, timser.length=max(index(timser.irreg)), int.method="linear", appliedFilter="lowpass", k=5)
##' plot(timser.irreg,type="l",col="black")
##' lines(timser.irregSmooth,col="limegreen",lwd=3)
##' @author
##' @export



InterpolationMethod <- function(X, fc, dt, timser.length, int.method=c("linear","nearest"), appliedFilter=c("gauss","runmean","lowpass"), k=5)
{
	if (!is.zoo(X)){stop("time series must be zoo-object")}
	#
	tx <- index(X)
	x <- coredata(X)
	int.out <- seq(from=dt, to=max(tx), by=dt)
	Xsmooth <- rep.int(NA, length(seq(from=dt, to=timser.length, by=dt)))
	#
	if(int.method=="linear"){ X.int <- zoo(approx(tx, x, int.out, method="linear")$y, order.by=int.out) }
	if(int.method=="nearest"){ X.int <- approx.nearest(tx, x, int.out) }	
	#
		if(appliedFilter=="gauss"){
			fw <- gauss.weights(fc)
			filt.pos <- seq(from=-(round(length(fw)/(2*dt))*dt), to=round(length(fw)/(2*dt))*dt, by=dt) + 1 + (length(fw)-1)/2
				if(max(filt.pos) > length(fw)){ filt.pos <- seq(from=-(round(length(fw)/(2*dt))-1)*dt, to=(round(length(fw)/(2*dt))-1)*dt, by=dt)+ 1 + (length(fw)-1)/2 }
			#
			Xsmooth[1:length(int.out)] <- filter(X.int, fw[filt.pos]/sum(fw[filt.pos]))
		}
		if(appliedFilter=="runmean"){
			fw <- runmean.weights(fc)
			filt.pos <- seq(from=-(round(length(fw)/(2*dt))*dt), to=round(length(fw)/(2*dt))*dt, by=dt) + 1 + (length(fw)-1)/2
				if(max(filt.pos) > 1/(2*fc)){ filt.pos <- seq(from=-(round(length(fw)/(2*dt))-1)*dt, to=(round(length(fw)/(2*dt))-1)*dt, by=dt) + 1 + (length(fw)-1)/2 }
			#
			Xsmooth[1:length(int.out)] <- filter(X.int, fw[filt.pos]/sum(fw[filt.pos]))
		}
		if(appliedFilter=="lowpass"){
			fw <- lowpass.weights(fc)
			filt.pos <- seq(from=-(round(length(fw)/(2*dt))*dt), to=round(length(fw)/(2*dt))*dt, by=dt) + 1 + (length(fw)-1)/2
				if(max(filt.pos) > length(fw)){ filt.pos <- seq(from=-(round(length(fw)/(2*dt))-1)*dt, to=(round(length(fw)/(2*dt))-1)*dt, by=dt)+ 1 + (length(fw)-1)/2 }
			#
			Xsmooth[1:length(int.out)] <- filter(X.int, fw[filt.pos]/sum(fw[filt.pos]))
		}
	#
	return(zoo(Xsmooth, order.by=seq(from=dt, to=timser.length, by=dt)))	
}


