#' Make tables with measurement only in the middle of interval.
#' @param tables - List of tables that were created by \emph{calaculate_tables} function.

#' @export

collapse.to.middle <- function(tables){
  tables.m <- list()
  
  for(i in 1:length(tables)){
    tV <- tables[i][[1]][1]$OD
    tC <- tables[i][[1]][2]$Conc
    tube <- names(tables)[i]
    tV.m <- data.frame(time = c(), tube = c(), value = c(), pump = c())    #data frames for middle point concentration and OD
    tC.m <- data.frame(time = c(), Concentration = c())
    for( j in 1:(length(tC$time)-1) ){
      tV.temp <- tV[tV$time >= tC$time[j] & tV$time < tC$time[j+1],] #take interval between dilutions
      if( length(tV.temp$time) < 2 ){
        next
      }
      OD.temp <- (tV.temp$value[1] + tail(tV.temp$value,1)) / 2           #assume that interval should be linear and take middle of the line
      Conc.temp <- tC$Concentration[j]
      tV.m <- rbind(tV.m, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/2 ), tube = tube, value = OD.temp, pump = tV.temp$pump[1]))
      tC.m <- rbind(tC.m, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/2 ), Concentration = Conc.temp))
    }
    
    #add last point
    tV.temp <- tV[tV$time >= tC$time[length(tC$time)],]
    if( length(tV.temp$time) > 0 ){
      OD.temp <- (tV.temp$value[1] + tail(tV.temp$value,1)) / 2
      Conc.temp <- tC$Concentration[length(tC)]
      tV.m <- rbind(tV.m, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/2 ), tube = tube, value = OD.temp, pump = tV.temp$pump[1]))
      tC.m <- rbind(tC.m, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/2 ), Concentration = Conc.temp))
    }
    
    table.m.names <- names(tables.m)
    tables.m <- rlist::list.append(tables.m, x = list(OD = tV.m,Conc = tC.m))
    names(tables.m) <- c(table.m.names, tube)
  }
  
  return(tables.m)
}