library("RColorBrewer")

#Setup directory list
resultsDirectory="/home/kathryn/Documents/Chaste/testoutput/"
zoneLengths = c(20,45,70)
springs = c(15,30,45,60)
resists = c(30,60,90,120)
parameterTest = c()
for(zl in zoneLengths){
for(sp in springs){
for(re in resists){
	parameterTest = c(parameterTest, paste("POSITIONAL_",zl,"_",sp,"_",re,sep=""))
}
}
}
dirs = parameterTest


#Setup for plotting
N = length(dirs)
par( mfrow=c(ceiling(sqrt(N)),ceiling(sqrt(N))), mar=c(1,1,1,1) )
includeDiffCells = FALSE
palette = brewer.pal(7, "Set1")
colorcode = c(palette[1], palette[5], palette[3], palette[2], palette[4], "darkblue","black")


#Processing and plotting loop
for( item in seq(1,N,by=1) ){

	phasedata<-read.table(paste(resultsDirectory, dirs[item], "/results_from_time_0/PhaseData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
	
	spl = strsplit(dirs[item],"_")
	z = spl[[1]][2] 
	s = spl[[1]][3] 
	r = spl[[1]][4] 
	header = paste("Resist: ", r, " Zone: ", z, " Spring: ", s)

	times = phasedata[,1]
	G1 = phasedata[,2] + phasedata[,6]
	S = phasedata[,3]
	G2 = phasedata[,4]
	M = phasedata[,5]
	MeiS = phasedata[,7]
	Diff = phasedata[,8]
	if(includeDiffCells){
	    Total = G1+S+G2+M+MeiS+Diff
	}else{
	    Total = G1+S+G2+M+MeiS
	}
	
	plot(times,   100*(G1/Total), type='l', lwd=2, main=header, xlab="Time (hph)", ylab="Percentage cells", ylim=c(0,100), col=colorcode[1])
	points(times, 100*(S/Total), type='l', lwd=2, col=colorcode[2])
	points(times, 100*(G2/Total), type='l', lwd=2, col=colorcode[3])
	points(times, 100*(M/Total), type='l', lwd=2, col=colorcode[4])
	points(times, 100*(MeiS/Total), type='l', lwd=2, col=colorcode[5])
	
	if(includeDiffCells){
	    points(times, 100*(Diff/Total), type='l', lwd=2, col=colorcode[6])
	    #legend("topright", legend=c("G1","S","G2","M","MeiS","Diff"), fill=colorcode)
	}else{
	    #legend("topright", legend=c("G1","S","G2","M","MeiS"), fill=colorcode)
	}
}	