#' Calculate number of events of replication
#' @param N.start - Number of cells at the beginning of the interval.
#' @param N.final - Number of cells at the end of the interval. 

#' @export

num.rep.events <- function(N.start, N.final){
    return(N.final - N.start)
}