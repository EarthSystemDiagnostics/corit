##' Simulate a random timeseries with a powerlaw spectrum
##'
##' Method: FFT white noise, rescale, FFT back, the result is scaled to variance 1
##' @title Simulate a random timeseries with a powerlaw spectrum
##' @param beta slope
##' @param N length of timeseries to be generated
##' @return vector containing the timeseries
##' @author Thomas Laepple
##' @export


SimPowerlaw <- function(beta, N)
{
  N2 <- (3^ceiling(log(N, base = 3)))
  df  <- 1 / N2
  f <- seq(from = df, to = 1/2, by = df)
  Filter <- sqrt(1/(f^beta))
  Filter <- c(max(Filter), Filter, rev(Filter))
  x   <- scale(rnorm(N2, 1))
  fx  <- fft(x)
  ffx <- fx * Filter
  result <- Re(fft(ffx, inverse = TRUE))[1:N]
  return(scale(result)[1:N])
}
