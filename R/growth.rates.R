#' Calculate list of table for growth rates (h^-1) at each interval with concentration on each interval, start and stop time for interval.
#' @param tables - List of tables that were created by \emph{calaculate_tables} function.
#' @param V - Volume of of media in the reactors.
#' @param MIC - Minimum Inhibitory Concentration.

#' @export

growth.rates <- function(tables, V = 20, MIC = 0.016){
  #V - volume of the reaction
  tables.GR <- list()
  for(i in 1:length(tables)){
    tV <- tables[i][[1]][1]$OD
    #tV$value <- tV$value - min(tV$value)
    tV$value <- log(tV$value)    #convert to natural logarithm of OD
    #tV$time <- tV$time/3600 #convert to hours
    tC <- tables[i][[1]][2]$Conc
    tC$Concentration <- tC$Concentration / MIC   #Concentration will be in the xMIC
    tube <- names(tables)[i]
    
    tG <- data.frame(time = c(), growth.rate = c(), Concentration = c(), t.start = c(), t.stop = c()) #Growth rates will be calculated for the middle of the interval
    
    for( j in 1:(length(tC$time)-1) ){
      tV.temp <- tV[tV$time >= tC$time[j] & tV$time < tC$time[j+1],]  #take interval betwee dilutions
      if( length(tV.temp$time) < 2 ){  #if interval has no values skip it
        #tG <- rbind(tG, data.frame( time = ( (tC$time[j] + tC$time[j+1])/7200 ), growth.rate = 0, Concentration = tail(tG$Concentration,1), t.start = tail(tV$time[tV$time < tC$time[j]],1)/3600, t.stop = head(tV$time[tV$time > tC$time[j+1]],1)/3600))
        next
      }
      
      fit <- lm(value ~ time, data = tV.temp)
      Conc.temp <- tC$Concentration[j]
      if( fit$coefficients[2] >= 0 ){
        tG <- rbind(tG, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/7200 ), growth.rate = fit$coefficients[2] * 3600, Concentration = Conc.temp, t.start = tV.temp$time[1]/3600, t.stop = tail(tV.temp$time,1)/3600))
      }else{
        tG <- rbind(tG, data.frame( time = ( (tV.temp$time[1] + tail(tV.temp$time,1))/7200 ), growth.rate = 0, Concentration = Conc.temp, t.start = tV.temp$time[1]/3600, t.stop = tail(tV.temp$time,1)/3600))
      }
    }
    
    tables.GR.names <- names(tables.GR)
    tables.GR <- rlist::list.append(tables.GR, x = tG)
    names(tables.GR) <- c(tables.GR.names, tube)
  }
  
  return(tables.GR)
}