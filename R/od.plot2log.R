#' Convert plot save format to log format.
#' @param plotfile - file with plot save. default = OD_plot.txt
#' @param logfile - path where log format should be saved. default = OD_log.txt

#' @export

od.plot2log <- function(plotfile = "OD0.txt", logfile = "OD_log.txt"){
	od.plot = read.table(plotfile, sep = '\t', header = T)
	od.log = data.frame(time=c(), od=c())
	
	#iterate throw column pairs date/OD
	for(i in 1:(dim(od.plot)[2]/2)){
	  tube.pump = strsplit(names(od.plot)[i*2-1], "[.]")
	  
	  if(tube.pump[[1]][2] == 'P1' ){
	    tube.pump[[1]][2] = 'Pump1'
	  }else{
	    tube.pump[[1]][2] = 'Pump2'
	  }
	  
	  od.plot.temp = od.plot[,c(2*i-1, 2*i)]
	  od.plot.temp = od.plot.temp[!is.na(od.plot.temp[,2]),]
	  
	  od.log.temp = data.frame(
	    time = od.plot.temp[,1], 
	    od = paste(tube.pump[[1]][1], od.plot.temp[,2], tube.pump[[1]][2], sep = ',')
	  )
	  
	  od.log = rbind(od.log, od.log.temp)
	}
	
	od.log = od.log[order(as.POSIXct(od.log$time)),]
	write.table(od.log, logfile, sep = '\t', col.names = F, row.names = F, quote = F)
}