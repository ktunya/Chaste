library("RColorBrewer")

resultsDirectory="/Users/kathryn/Code/HasteRemote/testoutput/"
dir = "Test50HSCI"

posdata<-read.table(paste(resultsDirectory, dir, "/results_from_time_0/PositionData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");

times = unique(posdata[,1])
timeDiscretisation = times[seq(1,length(times),by=20)]
spatialDiscretisation = seq(-15,250,by=50)
lenSp = length(spatialDiscretisation)

for(s in seq(1,lenSp-1,by=1) ){

	S1=spatialDiscretisation[s]
	S2=spatialDiscretisation[s+1]

	results = rep(0,length(timeDiscretisation)-1)

	for(i in seq(1,length(timeDiscretisation)-1,by=1)){
	
		T1=timeDiscretisation[i]
		T2=timeDiscretisation[i+1]
	
		selectedCellsT1 = posdata[(posdata[,1]==T1 & posdata[,3]>S1 & posdata[,3]<S2),]
		selectedCellsT2 = posdata[(posdata[,1]==T2 & posdata[,3]>S1 & posdata[,3]<S2),]
		cellIds = selectedCellsT1[selectedCellsT1[,2]%in%selectedCellsT2[,2],2];
	
		P1 = selectedCellsT1[selectedCellsT1[,2]%in%cellIds,3:5]
		P2 = selectedCellsT2[selectedCellsT2[,2]%in%cellIds,3:5]
	
		d = dim((P2-P1))
		meanComponents = colSums((P2-P1))/(d[1]*(T2-T1))
		results[i] = meanComponents[1]
	}

	if(s==1){
		plot(timeDiscretisation[2:length(timeDiscretisation)], results, col=rgb(1,0,0,1), type='l', lwd=2, xlab="Time (hph)", ylab="Mean cell speed",ylim=c(0,30))
	}else{
		points(timeDiscretisation[2:length(timeDiscretisation)],results, col=rgb(1-s/lenSp,0,s/lenSp,1), type='l', lwd=2)
	}
}


