##'Averages y into bins according to the positon of a in the breaks
##'Either give N=number of breaks, or N+1 breaks
##' Breaks are defined as x>breaks[i], and x<=breaks[i+1]
##' if fill=T, fill empty bins using linear interpolation from the neighbours to the center of the bin
##' could be considerably speeded up by using cut ?cut
##' @title Average a vector into bins
##' @param x vector of x values
##' @param y vector of y values, same length as x
##' @param N Number of breaks (or NULL if breaks are supplied)
##' @param breaks vector of breaks (optional, instead if N)
##' @param bFill if TRUE  fill empty bins using linear interpolation from the neighbours to the center of the bin
##' @return list(breaks,centers,avg,nobs)
##' Returns the breaks, centers, the averaged values and nobs, the number of observations averages
##' @author Thomas Laepple
##' @export


AvgToBin<-function(x,y,N=NULL,breaks=pretty(x,N),bFill=FALSE)
{
    NBIN <- length(breaks) - 1
    centers <- (breaks[1:NBIN] + breaks[2:(NBIN + 1)])/2
    avg<-rep(NA,length(breaks)-1)
    nobs<-rep(NA,length(breaks)-1)
    for (i in 1:(length(breaks)-1)) {
        selection<-y[which((x>breaks[i])&(x<=breaks[i+1]))]
        avg[i]<-mean(na.omit(selection))
        nobs[i]<-sum(!is.na(selection))
        }


  if ((sum(is.na(avg))>0)&(bFill))
  {
      yInt<-approx(x,y,centers)$y
      missing<-is.na(avg)
      avg[missing]<-yInt[missing]
  }
   return(list(breaks=breaks,centers=centers,avg=avg,nobs=nobs))
}
