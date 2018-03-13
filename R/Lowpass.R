##' Derive the (smoothed) least square lowpass, given the cutoff frequency omega.c and the length of the filter n
##' based on Bloomfield 1976
##' 
##' @title calculate weights for lowpass filter
##' @param omega.c cutoff frequency
##' @param n length of the filter, has to be odd
##' @param sample sampling rate of the timeseries on which the filter will be applied (1/deltat)
##' @param convergence TRUE: smoothed least square lowpass; FALSE = unsmoothed
##' @return vector of filter weights
##' @author Thomas Laepple
##' @export


Lowpass <- function(omega.c,
                    n = 9,
                    sample = 1,
                    convergence = TRUE){
  if ((n %% 2) == 0)
    stop("N must be odd, this function calculates only symmetrical = phase preserving filters")
  
  omega.c <- omega.c / sample
  if (omega.c >= 0.5)
    stop("frequency higher or equal then Nyquist")
  #calculate least square solution for ideal lowpass
  n.side <- (n - 1) / 2
  g.u <- numeric(length = n)
  
  u <- 1:n.side
  g.u[u + n.side + 1] <- sin(u * omega.c * 2 * pi) / (pi * u)
  g.u[n.side + 1] <- omega.c * 2
  #Now flip over
  g.u[1:n.side] <- g.u[n:(n.side + 2)]
  
  if (convergence) {
    #Multiply by convergence factors to reduce the ripples
    delta <- 4 * pi / n
    u <- c((-1 * n.side):-1, 1:n.side)
      g.u[u + n.side + 1] <- g.u[u + n.side + 1] * sin(u * delta / 2) / (u * delta / 2)
  }
  
  return(g.u)
}
