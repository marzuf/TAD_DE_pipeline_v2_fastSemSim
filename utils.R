printAndLog <- function(txt, logFile){
  cat(txt)
  cat(txt, file = logFile, append=T)
}

plot_multiDens <- function(size_list, plotTit="", legTxt=NULL, legPos="topright", my_ylab="density", my_xlab="") {
  
  dens <- lapply(size_list, function(x) density(na.omit(x)))
  names(dens) <- names(size_list)
  
  lengthDens <- unlist(lapply(size_list, function(x) length(na.omit(x))))
  
  plot(NA, xlim=range(sapply(dens, "[", "x")), ylim=range(sapply(dens, "[", "y")), 
       main=plotTit, xlab=my_xlab, ylab=my_ylab)
  foo <- mapply(lines, dens, col=1:length(dens))
  if(is.null(legTxt)){
    # legTxt <- names(dens)
    legTxt <- paste0(names(dens), " (n=", lengthDens, ")")
  }
  legend(legPos, legend=legTxt, fill=1:length(dens), bty='n')
}

######################################################################################################################################################################################################
###################################################################################################################################################################################################### densplot (function)
######################################################################################################################################################################################################

densplot <- function(x,y, pch=19, cex=1, ...){
	df <- data.frame(x,y)
	d <- densCols(x,y, colramp=colorRampPalette(c("black", "white")))
	df$dens <- col2rgb(d)[1,] + 1L
	cols <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F","#FCFF00", "#FF9400", "#FF3100"))(256)
	df$col <- cols[df$dens]
	df <- df[order(df$dens),]
	plot(df$x,df$y, pch=pch, col=df$col, ...)
}

