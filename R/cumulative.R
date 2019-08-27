#' Cumulative values of OD, Concentration, number of generations and replication events.
#' @param tables - List of tables that were created by \emph{calaculate_tables} function.
#' @param KOD - coefficient of the OD to number of cell conversion.
#' @param V - volume of media in the reactors.

#' @export

cumulative <- function(tables, KOD = 8 * 10^8, V = 20){
  tables.cum <- list()
  
  for(i in 1:length(tables)){
    tV <- tables[i][[1]][1]$OD
    tC <- tables[i][[1]][2]$Conc
    tube <- names(tables)[i]
    OD.cum <- tV$value[1]
    Conc.cum <- tC$Concentration[1]
    
    tV.cum <- data.frame(time = tV$time[1], tube = tube, value = OD.cum, pump = tV$pump[1], Gen.cum = 0, N.rep.events = 0)    #data frames for middle point concentration and OD + number of generations from the time 0
    tC.cum <- data.frame(time = tV$time[1], Concentration = Conc.cum)
    
    tV.temp <- tV[tV$time < tC$time[1],]      #first interval
    if(length(tV.temp$time) > 1){
      fit <- lm(value ~ time, data = tV.temp)
      OD.stop <- predict(fit,data.frame(time = tail(tV.temp$time,1)))
      OD.start <- predict(fit,data.frame(time = tV.temp$time[1]))
      tV.cum$Gen.cum[1] <- log(OD.stop/OD.start,base = 2)
      tV.cum$N.rep.events[1] <- num.rep.events(OD.start * V * KOD, OD.stop * V * KOD )
    }
    
    for( j in 1:(length(tC$time)-1) ){
      tV.temp <- tV[tV$time >= tC$time[j] & tV$time < tC$time[j+1],] #take interval between dilutions
      
      if( length(tV.temp$time) <= 1 ){  #if interval has no values skip it
        next
      }
      
      fit <- lm(value ~ time, data = tV.temp)
      #OD.temp <- (tail(tV.temp$value,1) - tV.temp$value[1])
      OD.temp <- fit$coefficients[2] * (tail(tV.temp$time,1) - tV.temp$time[1])
      
      if( OD.temp < 0 ){
        OD.temp <-  0
        Gen.cum.temp <- tail(tV.cum$Gen.cum,1)
        N.rep.events.temp <- tail(tV.cum$N.rep.events,1)
      }else{
        OD.stop <- predict(fit,data.frame(time = tail(tV.temp$time,1)))
        OD.start <- predict(fit,data.frame(time = tV.temp$time[1]))
        if( OD.stop/OD.start > 1 ){   # culture is growing is or not. If not, then just take last generations calculations
          Gen.cum.temp <- tail(tV.cum$Gen.cum,1) + log(OD.stop/OD.start,base = 2)
          N.rep.events.temp <- tail(tV.cum$N.rep.events,1) + num.rep.events(OD.start * V * KOD, OD.stop * V * KOD )
        }else{
          Gen.cum.temp <- tail(tV.cum$Gen.cum,1)
          N.rep.events.temp <- tail(tV.cum$N.rep.events,1)
        }
      }
      
      OD.cum <- OD.cum + OD.temp
      Conc.cum <- Conc.cum + tC$Concentration[j]
      tV.cum <- rbind(tV.cum, data.frame( time = ((tail(tV.temp$time,1) + tV.temp$time[1])/2), tube = tube, value = OD.cum, pump = tV.temp$pump[1], Gen.cum = Gen.cum.temp, N.rep.events = N.rep.events.temp))
      tC.cum <- rbind(tC.cum, data.frame( time = ((tail(tV.temp$time,1) + tV.temp$time[1])/2), Concentration = Conc.cum))
    }
    
    #add last point
    tV.temp <- tV[tV$time >= tC$time[length(tC$time)],]
    if( length(tV.temp$time) > 1 ){
      fit <- lm(value ~ time, data = tV.temp)
      if( fit$coefficients[2] < 0 ){
        Gen.cum.temp <- tail(tV.cum$Gen.cum,1)
        N.rep.events.temp <- tail(tV.cum$N.rep.events,1)
      }else{
        OD.stop <- predict(fit,data.frame(time = tail(tV.temp$time,1)))
        OD.start <- predict(fit,data.frame(time = tV.temp$time[1]))
        Gen.cum.temp <- tail(tV.cum$Gen.cum,1) + log(OD.stop/OD.start,base = 2)
        N.rep.events.temp <- tail(tV.cum$N.rep.events,1) + num.rep.events(OD.start * V * KOD, OD.stop * V * KOD )
      }
      OD.temp <- tail(tV.temp$value,1) - tV.temp$value[1]
      Conc.cum <- Conc.cum + tC$Concentration[length(tC$time)]
      tV.cum <- rbind(tV.cum, data.frame( time = ((tail(tV.temp$time,1) + tV.temp$time[1])/2), tube = tube, value = OD.cum, pump = tV.temp$pump[1], Gen.cum = Gen.cum.temp, N.rep.events = N.rep.events.temp))
      tC.cum <- rbind(tC.cum, data.frame( time = ((tail(tV.temp$time,1) + tV.temp$time[1])/2), Concentration = Conc.cum ))
    }
    
    table.cum.names <- names(tables.cum)
    tables.cum <- rlist::list.append(tables.cum, x = list(OD = tV.cum,Conc = tC.cum))
    names(tables.cum) <- c(table.cum.names, tube)
  }
  
  return(tables.cum)
}