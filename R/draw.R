#' Draw OD/time and Concentration/time curves.
#' @param tables - List of tables that were created by \emph{calaculate_tables} function.
#' @param c.units - Text variable that specifies concentration units.
#' @param MIC - Minimum Inhibitory Concentration.
#' @param log2 - Should the Concentration axes transformed with log2 transformation? Default value is TRUE (log2 coordinates).
#' @param p.size - Size of points on the OD plot.
#' @param add.model - If mutation dynamics modelling was made with the \emph{mutability} function TRUE value will add this model to OD plot. Default value is FALSE.
#' @param model.table - Table with the model dynamics produced by \emph{mutability} function. Could be omitted if \strong{add.model == FALSE}.
#' @param m.prob - Mutation probability. Could be omitted if \strong{add.model == FALSE}. Calculated inside the \emph{calaculate_tables} function.
#' @param params - information about dilutions and bottle changes geneated by \emph{get_draw_parameters}
#' @param dist - distance between lines indicated sample collections or bottle changes and text
#' @param ymeta - height of bottle change and sample collection lables

#' @export

#' @import ggplot2
#' @import gtable
#' @import grid

draw <- function(tables, c.units="uM", MIC=0.016, log2=T, p.size = 0.5, add.model = F, model.table = data.frame(), m.prob = "", params = list(), dist = 2, ymeta = 0.3){
  #if log2==T, then concentration is shown as log2 of xMIC. Else - in linear coordinates in c.units
  #tables - list of lists with tables of concentration and OD
  #p.size - size of points on the OD plot
  if( length(params) != 0 ){
    params[['bcs']]$time = params[['bcs']]$time / 3600
    params[['dils']]$time = params[['dils']]$time / 3600
  }

  for(i in 1:length(tables)){
    tV <- tables[i][[1]][1]$OD
    tC <- tables[i][[1]][2]$Conc
    maxCS <- max(tC$Concentration)
    tC$time <- tC$time/3600   #converts time to hours
    tV$time <- tV$time/3600
    tube <- names(tables)[i]
    if( log2 == F){
      tplotC <- ggplot(tC,aes(time, Concentration)) +
        geom_line() +
        ggtitle("Antibiotic concentration\n") +
        xlim(min(tV$time,tC$time)-4,max(tV$time,tC$time)) +
        labs(y = paste("Concentration",c.units),x = "Time, h") +
        theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9), plot.title = element_text(face="bold")) +
        ylim(0,maxCS) +
        theme_bw()
    }else{
      tplotC <- ggplot(tC[tC$Concentration > MIC/10,],aes(time, Concentration/MIC)) +
        geom_line() +
        ggtitle("Antibiotic concentration\n") +
        xlim(min(tV$time,tC$time)-4,max(tV$time,tC$time)) +
        labs(y = bquote('xMIC ('~log[2]~' scale)'),x = "Time, h") +
        theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9), plot.title = element_text(face="bold")) +
        ylim(0,maxCS/MIC) +
        scale_y_continuous(trans='log2',breaks = geomSeries(base=2, max=maxCS/MIC)) +
        theme_bw()
    }
    tC$panel <- "Concentration"

    tV$panel <- "OD"
    tplotV <- ggplot(tV,aes(time, value)) +
      geom_line() +
      geom_point(aes(color = pump), size=p.size) +
      ggtitle(paste0("OD ",tube,"\n")) +
      xlim(min(tV$time,tC$time)-4,max(tV$time,tC$time)) +
      labs(y = "OD",x = "Time, h") +
      theme(axis.title.x = element_text(size=9),axis.title.y = element_text(size=9), plot.title = element_text(face="bold")) +
      theme_bw() +
      scale_color_manual(values = c("green2","red2"))

    if( add.model == T){
      model.table$time.mut = model.table$time.mut/3600
      tplotV <- tplotV +
        geom_line(data = model.table, aes(time.mut, value.mut),colour = "blue") +
        geom_label(aes(x = 50, y = 0.1, label = paste("P(mutation) =", formatC(m.prob, format = "e", digits = 2))))
    }

    if( length(params) != 0 ){
      for( b in 1:dim(params[['bcs']])[1] ){
        tplotV <- tplotV +
          geom_vline(xintercept = params[['bcs']]$time[b], colour = "blue") +
          annotate("text", x=params[['bcs']]$time[b] - dist, label=paste0("M1 = ", round(params[['bcs']]$Cp1[b]/MIC, 0), "x, M2 = ", round(params[['bcs']]$Cp2[b]/MIC, 0), "x"), y=ymeta, colour="blue", angle=90, size=4)

        tplotC <- tplotC +
          geom_vline(xintercept = params[['bcs']]$time[b], colour = "blue")
      }

      temp_dils =  params[['dils']][params[['dils']]$tube == tube,]

      if(length(temp_dils) > 0){
      	for( d in 1:dim(temp_dils)[1] ){
      	  tplotV <- tplotV +
      	    geom_vline(xintercept = temp_dils$time[d], colour = "red") +
      	    annotate("text", x=temp_dils$time[d] + dist, label=paste0("Sample ", round(temp_dils$time[d], 0), "hrs"), y=ymeta, colour="red", angle=90, size=4)
	
      	  tplotC <- tplotC +
      	    geom_vline(xintercept = temp_dils$time[d], colour = "red")
      	}
      }

    }

    gC <- ggplotGrob(tplotC)
    gV <- ggplotGrob(tplotV)
    for(j in 1:(ncol(gV) - ncol(gC))){
      gC <- gtable::gtable_add_cols(gC, unit(0,"mm"))
    }
    g <- rbind(gV, gC, size="first") # stack the two plots
    g$widths <- unit.pmax(gC$widths, gV$widths)
    g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
    grid.newpage()
    pdf(paste0(tube,".pdf"))
    grid.draw(g)
    dev.off()
  }
}
