library("RColorBrewer")

resultsDirectory="/Users/kathryn/Code/HasteRemote/testoutput/"
dir = "Test50HSCI"

phasedata<-read.table(paste(resultsDirectory, dir, "/results_from_time_0/MitosisData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");

times = unique(phasedata[,1])
percentRange = c(80,99)
range = times[floor((percentRange[1]/100)*length(times)) : floor((percentRange[2]/100)*length(times))]

total = rep(0,50)
mitosis = rep(0,50)

for(t in range){

	selectedCells = phasedata[phasedata[,1]==t,]
	selectedCells = selectedCells[2:length(selectedCells)]
	nRows = length(selectedCells)/2

	total = total + as.vector(selectedCells[ 1:nRows ])
	mitosis = mitosis + as.vector(selectedCells[ (nRows+1):(2*nRows) ])
}

MI = 100*(mitosis/total)

#c=rgb((t/last),0,(1-t/last),1)
#if(t==times[1]){
plot(seq(1,nRows,by=1), MI, col="black", type='l', lwd=2, xlab="Cell row", ylab="Mitotic Index (%)", ylim=c(0,5))
#}else{
#	points(seq(1,nRows,by=1), MI, col=c, type='l', lwd=2)
#}