##' @title Fit a power-law to the spectrum
##' @param spec
##' @param freq.start  vector containing the start frequencies of the fitting interval(s)
##' @param freq.end vector containing the end frequencies of the fitting interval(s)
##' @param bDebug (TRUE) plot diagnostics
##' @param breaks  vector of breakpoints to which the spectra is binned (optional)
##' @param indexRemove bins that are removed (e.g. containing the annual and semiannual cycle)
##' @param df.log resolution of the bins (if breaks are not provided)
##' @param i.fStart index of first (lowest) frequency to be used
##' @return list(slope=slope,slopesd=slopesd,spec=saveSpec,freq=binFreq,intercept=intercept)
##' @author Thomas Laepple
##' @export


SlopeFit<-function(spec,freq.start=NULL,freq.end=NULL,bDebug=TRUE,breaks=NULL,indexRemove=NULL,df.log=0.05,i.fStart=4)
{
   if (is.null(breaks)) breaks=seq(from=log(spec$freq[i.fStart]),to=log(tail(spec$freq, 1L)),by=df.log)
   binSpec<-AvgToBin(spec$freq,log(spec$spec),N=nBins,bFill=T,breaks=exp(breaks))
#Remove missing values (if the breaks extended over the frequency range
    index<-!is.na(binSpec$avg)
    binSpec$avg<-exp(binSpec$avg)
    binSpec$avg[indexRemove]<-0
    binSpec$avg<-binSpec$avg[index]

    binSpec$centers<-binSpec$centers[index]

    saveSpec<-binSpec$avg



    binFreq<-binSpec$centers #Get back of the
    binSpec$centers=log(binSpec$centers)




    if (is.null(freq.start)) freq.start<-min(binFreq)
    if (is.null(freq.end)) freq.end  <-max(binFreq)
    if (length(freq.start) != length(freq.end)) stop("The same number of start and stop frequencies are needed")
    slope<-slopesd<-NA



if (bDebug)
{
    plot(spec$freq,spec$spec,log="xy",type="l",xlab="f (1/yr)",ylab="PSD")
    lines(binFreq,binSpec$avg,col="green",lwd=2)
}

#Remove the annual cycle and start and freq.end from the spectra


slope<-slopesd<-intercept<-NA
for (i in 1:length(freq.start))
{
    index<-rep(TRUE,length(binFreq))

    i.start<-ClosestElement(binFreq,freq.start[i]) #Remove the low frequency part
    i.end<-ClosestElement(binFreq,freq.end[i]) #Remove the high frequency part

    index[binSpec$avg==0]<-FALSE
    index[1:(i.start-1)]<-FALSE
    index[(i.end+1):length(index)]<-FALSE


    if (bDebug)
    {
        abline(v=binFreq[c(i.start,i.end)])
    }

#Linear fit...
    model<- lm(log(binSpec$avg[index])~binSpec$centers[index])

    if (bDebug)   lines(binFreq[index],exp(model$fitted.values),lwd=1,col="red")
    slopesd[i]<-summary(model)$coeff[2,2]
    slope[i]<-model$coeff[2]
     intercept[i]<-model$coeff[1]
}
return(list(slope=slope,slopesd=slopesd,spec=saveSpec,freq=binFreq,intercept=intercept))
}
