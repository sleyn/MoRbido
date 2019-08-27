#' Create a list tables of OD and Concentration changes through time from the log files.
#' @param dil_check - Were any manual dilutions introdused? Default is TRUE. If TRUE - need to provide dil_file
#' @param dil_file - File that specifies dilution time, ratio (volume taken to total volume), concentration of the drug in the added media and the reactor.
#' @param bc_check - Was concentration of a drug changed during the run. Default is TRUE.
#' @param bc_file - File that specifies bottle change time, concentration in the Pump 1 bottle, Pump 2 bottle.
#' @param C0 - Concentration in reactors at the start of the run.
#' @param CS1 - Concentration in Pump 1 bottle at the start of the run.
#' @param CS2 - Concentration in Pump 2 bottle at the start of the run.
#' @param V - Volume of of media in the reactors.
#' @param dV - The delution volume.
#' @param od_file - Log file of OD measurements.
#' @param pump_file - Log file of pumps activation.
#' @param tubes_selected - Vector with numbers of reactors activated in the current run.
#' @param zero - zero levels of each tube. The levels should be in the increasing order.
#' @param dV.change - file with time of change and new dV

#' @export

#' @import rlist

calculate_tables <- function(dil_check = T, bc_check = T, C0 = 0, CS1 = 0, CS2 = 10, V =20, dV = 4, od_file="6-21-17-OD.txt",pump_file="6-21-17-pumps.txt",dil_file="dilutions.txt",bc_file = "bottle_change.txt", tubes_selected = 1:6, zero = c(Tube1 = 0,Tube2 = 0,Tube3 = 0,Tube4 = 0,Tube5 = 0,Tube6 = 0), dV.change = c()){
  #C0 antibiotic concentration in tubes on start (ug/ml)
  #CS1 antibiotic concentration in media 1 (ug/ml)
  #CS2 antibiotic concentration in media 2 (ug/ml)
  #V reaction volume (ml)
  #dV dilution volume (ml)
  #dil_check True if dilutions were present
  #bc_check True if botle changes were present
  
  table.t <- read.table(od_file,sep="\t")      #Read OD table
  colnames(table.t) <- c("time","tube","value","pump")
  table.t$time <- as.POSIXlt(table.t$time)
  time0 <- table.t$time[1]
  table.t$time <- as.numeric(table.t$time - time0, units = "secs")    #convers time to seconds passed from the first OD measurement
  
  for( tube.name in names(zero)){			#introduce zero values
  	table.t[table.t$tube == tube.name,]$value = table.t[table.t$tube == tube.name,]$value - zero[tube.name]
  }
  
  table.t$value[table.t$value <= 0] = 0.01
  #table.t <- table.t[table.t$value > 0,]
  
  table.p <- read.table(pump_file,sep="\t")     #Pump work table
  colnames(table.p) <- c("time","tube","pump")
  table.p$time <- as.POSIXlt(table.p$time)
  table.p$time <- as.numeric(table.p$time - time0, units = "secs")    #convers time to minutes passed from the first OD measurement
  
  if( length(dV.change) > 0 ){
    names(dV.change) = c('time', 'dV', 'tube')
    dV.change$time = as.POSIXlt(dV.change$time)
    dV.change$time =  as.numeric(dV.change$time - time0, units = "secs")
  }
  
  tables <- list() #list where will be stored tV and tC tables for tubes
  
  if( dil_check == T ){
    dilutions <- read.table(dil_file,sep="\t") #load dilution table
    colnames(dilutions) <- c("time","dilution","conc","tube")
    dilutions$time <- as.POSIXlt(dilutions$time)
    dilutions$time <- as.numeric(dilutions$time - time0, units = "secs")
  }
  
  if( bc_check == T){
    change.bottles <- read.table(bc_file,sep="\t") #load bottle change
    colnames(change.bottles) <- c("time","Cp1","Cp2")
    change.bottles$time <- as.POSIXlt(change.bottles$time)
    change.bottles$time <- as.numeric(change.bottles$time - time0, units = "secs")
  }
  
  for(i in tubes_selected){
    tube = paste0("Tube",i)
    
    tp <- table.p[table.p$tube == tube,]
    
    if( dil_check == T){
      d <- dilutions[dilutions$tube == tube,]
      dc <- 1 #dilutions count
    }
    
    if( bc_check == T ){
      bc <- 1 #bottle change count
    }
    
    Cn <- C0
    CS11 <- CS1
    CS21 <- CS2
    tC <- data.frame(c(),c(),colnames(c("time","Concentration")))
    ttemp <- data.frame(time = 0,Concentration = Cn)       #zero time point for concentration
    tC <- rbind(tC,ttemp)
    tV <- table.t[table.t$tube == tube,]
    
    if( length(dV.change) > 0 ){  # make change of dV table specific to the current tube
      dV.c.t = dV.change[dV.change$tube == tube, ]
      dV.c = 1      # counts of dV changes
    }
    
    if( dil_check == T && bc_check == T){   #BOTH Dilutions and Bottle changes were present
      if(dim(d)[1] > 0 ){
        for(j in 1:dim(d)[1]){          #clear jumps at dilution times
          tV <- tV[-which(tV$time > d$time[j] & tV$time < d$time[j] + 3600),]
        }
      }
      
      for(j in 1:dim(tp)[1]){ 
        if(length(d$time) >= dc){             #introduce dilutions
          if(tp$time[j] >= d$time[dc]){
            Cn <- Cn * (1 - d$dilution[dc]) + d$conc[dc] * d$dilution[dc]
            ttemp <- data.frame(as.integer((tp$time[j-1]+tp$time[j])/2),Cn)
            colnames(ttemp) <- c("time","Concentration")
            tC <- rbind(tC,ttemp)
            dc <- dc + 1
          }
        }
        
        if(length(change.bottles$time) >= bc){        #introduce bottle change
          if(tp$time[j] >= change.bottles$time[bc]){
            CS11 <- change.bottles$Cp1[bc]
            CS21 <- change.bottles$Cp2[bc]
            bc <- bc + 1
          }
        }
        
        if( length(dV.change > 0)){     # introduce change of dV
          if( length(dV.c.t$time) >= dV.c ){
            if( tp$time[j] >= dV.c.t$time[dV.c] ){
              dV = dV.c.t$dV
              dV.c = dV.c + 1
            }
          }
        }
        
        if(tp$pump[j] == "Pump1"){            #main concentration changes
          Cn <- (Cn*(V - dV)+CS11*dV)/V
        }else{
          Cn <- (Cn*(V - dV)+CS21*dV)/V
        }
        tC <- rbind(tC,data.frame(time = tp$time[j], Concentration = Cn))
      }
    }else if(dil_check == T && bc_check == F){  #Only dilutions wer present
      for(j in 1:dim(d)[1]){          #clear jumps at dilution times
        tV <- tV[-which(tV$time > d$time[j] & tV$time < d$time[j] + 3600),]
      }
      
      for(j in 1:dim(tp)[1]){ 
        if(length(d$time) >= dc){             #introduce dilutions
          if(tp$time[j] >= d$time[dc]){
            Cn <- Cn * (1 - d$dilution[dc]) + d$conc[dc] * d$dilution[dc]       #dilution is (V-dV)/V
            ttemp <- data.frame(tp$time[j-1]+1,Cn)
            colnames(ttemp) <- c("time","Concentration")
            tC <- rbind(tC,ttemp)
            dc <- dc + 1
          }
        }
        
        if( length(dV.change > 0)){     # introduce change of dV
          if( length(dV.c.t$time) >= dV.c ){
            if( tp$time[j] >= dV.c.t$time[dV.c] ){
              dV = dV.c.t$dV
              dV.c = dV.c + 1
            }
          }
        }
        
        if(tp$pump[j] == "Pump1"){            #main concentration changes
          Cn <- (Cn*(V - dV)+CS11*dV)/V
        }else{
          Cn <- (Cn*(V - dV)+CS21*dV)/V
        }
        tC <- rbind(tC,data.frame(time=c(tp$time[j]),Concentration=c(Cn)))
      }
    }else if(dil_check == F && bc_check == T){    #Only bottle changes were present
      for(j in 1:dim(tp)[1]){ 
        if(length(change.bottles$time) >= bc){        #introduce bottle change
          if(tp$time[j] >= change.bottles$time[bc]){
            CS11 <- change.bottles$Cp1[bc]
            CS21 <- change.bottles$Cp2[bc]
            bc <- bc + 1
          }
        }
        
        if( length(dV.change > 0)){     # introduce change of dV
          if( length(dV.c.t$time) >= dV.c ){
            if( tp$time[j] >= dV.c.t$time[dV.c] ){
              dV = dV.c.t$dV
              dV.c = dV.c + 1
            }
          }
        }
        
        if(tp$pump[j] == "Pump1"){            #main concentration changes
          Cn <- (Cn*(V - dV)+CS11*dV)/V
        }else{
          Cn <- (Cn*(V - dV)+CS21*dV)/V
        }
        tC <- rbind(tC,data.frame(time=c(tp$time[j]),Concentration=c(Cn)))
      }
    }else{
      for(j in 1:dim(tp)[1]){
        
        if( length(dV.change > 0)){     # introduce change of dV
          if( length(dV.c.t$time) >= dV.c ){
            if( tp$time[j] >= dV.c.t$time[dV.c] ){
              dV = dV.c.t$dV
              dV.c = dV.c + 1
            }
          }
        }
        
        if(tp$pump[j] == "Pump1"){            #main concentration changes
          Cn <- (Cn*(V - dV)+CS11*dV)/V
        }else{
          Cn <- (Cn*(V - dV)+CS21*dV)/V
        }
        tC <- rbind(tC,data.frame(time=c(tp$time[j]),Concentration=c(Cn)))
      }
    }
    
    table.names <- names(tables)
    tables <- rlist::list.append(tables,x = list(OD = tV,Conc = tC))
    names(tables) <- c(table.names,tube)
  }
  
  return(tables)
}