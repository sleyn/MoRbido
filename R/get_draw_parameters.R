#' Create a list of auxillary parameters for draw function
#' @param dil_check - Were any manual dilutions introdused? Default is TRUE. If TRUE - need to provide dil_file
#' @param dil_file - File that specifies dilution time, ratio (volume taken to total volume), concentration of the drug in the added media and the reactor.
#' @param bc_check - Was concentration of a drug changed during the run. Default is TRUE.
#' @param bc_file - File that specifies bottle change time, concentration in the Pump 1 bottle, Pump 2 bottle.
#' @param CS1 - Concentration in Pump 1 bottle at the start of the run.
#' @param CS2 - Concentration in Pump 2 bottle at the start of the run.
#' @param od_file - Log file of OD measurements.

#' @export

get_draw_parameters <- function(dil_check = T, bc_check = T, CS1 = 0, CS2 = 10, od_file="6-21-17-OD.txt", dil_file="dilutions.txt", bc_file = "bottle_change.txt"){

  # read first line of OD file and take start time
  table_f <- file(od_file, 'r')
  time0 <- strsplit(readLines(od_file,n=1),'\t')[[1]][1]
  time0 <- as.POSIXlt(time0)
  close(table_f)
  change.bottles = data.frame(time = 0, Cp1 = CS1, Cp2 = CS2)
  dilutions = data.frame(time = c(), dilution = c(), conc = c(), tube = c())

  # read dilution file
  if( dil_check == T ){
    dilutions <- read.table(dil_file,sep="\t") #load dilution table
    colnames(dilutions) <- c("time","dilution","conc","tube")
    dilutions$time <- as.POSIXlt(dilutions$time)
    dilutions$time <- as.numeric(dilutions$time - time0, units = "secs")
  }

  #read bottle change file
  if( bc_check == T){
    change.bottles <- read.table(bc_file,sep="\t") #load bottle change
    colnames(change.bottles) <- c("time","Cp1","Cp2")
    change.bottles$time <- as.POSIXlt(change.bottles$time)
    change.bottles$time <- as.numeric(change.bottles$time - time0, units = "secs")
  }
  
  change.bottles = rbind(data.frame(time = 0, Cp1 = CS1, Cp2 = CS2),change.bottles)

  params = list(dils = dilutions, bcs = change.bottles)
  
  if( dil_check == F & bc_check == F ){
  	params = list()
  }

  return(params)
}
