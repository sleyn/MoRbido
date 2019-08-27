#' Calculate a series of numbers in a geometric progression
#' with a given "base" up to "max" value.
#' @param base - A base of the geometric progression
#' @param max - A number that geometric progression should not exceed.

#' @export

geomSeries <- function(base, max) {
  base^(0:floor(log(max, base)))
}