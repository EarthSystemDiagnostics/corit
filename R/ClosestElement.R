##' Get closest element of a vector
##'
##' @param xvector a vector of values
##' @param x the value to find the closest match to
##' @param type
##'
##' @return
##' @export
##'
##' @examples


ClosestElement <- function(xvector, x, type = "N")
{
  if (type == "N")
    return(which.min(abs(x - xvector)))

  if (min(diff(xvector)) < 0)
    stop("Vector must be monotonically increasing for the the methods M and L")
  if (type == "M")
    return(first(which(x <= xvector)))
  if (type == "L")
    return(last(which(x >= xvector)))
}
