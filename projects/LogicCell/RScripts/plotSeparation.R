resultsDirectory="/Users/kathryn/Code/HasteRemote/testoutput/"
dirs = c("DistalArmThresholdSpring15", "DistalArmThresholdHC05_0.1Adaptive")
colorcode = c("dodgerblue", "red3")

for( d in seq(1,length(dirs), by=1) ){ 

    posdata<-read.table(paste(resultsDirectory, dirs[d], "/results_from_time_0/PositionData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
    times = unique(posdata[,1])
    avSpacings = rep(0,length(times))
    sdSpacings = rep(0,length(times))
    
    for(t in seq(1,length(times),by=1)){
    
    	selectedCells = posdata[posdata[,1]==times[t],2:6]
    
    	#print(t)
        #print(sum(is.nan(selectedCells[,5])))
    	avAtThisTime = median(selectedCells[!is.nan(selectedCells[,5]),5])
    	#sdAtThisTime = sd(selectedCells[!is.nan(selectedCells[,5]),5])
    	avSpacings[t] = avAtThisTime;
    	sdSpacings[t] = sdAtThisTime;
    
    	#if(t==length(times)-2){
    	#	hist(selectedCells[!is.nan(selectedCells[,5]),5])
    	#}
    	
    }
    
    if(d==1){
        plot(times, avSpacings, type="l", lwd=2, xlab="Time (hph)", ylab="Mean neighbour separation (microns)",
             ylim=c(min(avSpacings-sdSpacings), max(avSpacings+sdSpacings)), col=colorcode[d] )
    }else{
    	points(times, avSpacings, type="l", lwd=2, xlab="Time (hph)", ylab="Mean neighbour separation (microns)",
             ylim=c(min(avSpacings-sdSpacings), max(avSpacings+sdSpacings)), col=colorcode[d] )
    }
    #points(times, avSpacings+sdSpacings, type="l", lty=3, col=colorcode[d])
    #points(times, avSpacings-sdSpacings, type="l", lty=3, col=colorcode[d])

}