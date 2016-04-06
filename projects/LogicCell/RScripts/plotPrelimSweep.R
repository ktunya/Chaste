library("RColorBrewer")
library("plotrix")

resultsDirectory="/home/kathryn/Documents/Pos/"

test = c("TestPosCI")
pos = c("POSITIONAL_100_0.2_35_")#, "POSITIONAL_100_0.2_50_", "POSITIONAL_100_0.2_70_")
posCI = c("CIPOSITIONAL_100_0.2_35_", "CIPOSITIONAL_100_0.2_50_", "CIPOSITIONAL_100_0.2_70_")
posSlow = c("SLOWPOSITIONAL_100_0.2_35_10_")
gen = c()
genCI = c()
genSlow = c()
dirs = posCI
repNums = seq(0,8,by=1);


#mitosis options / initialisation
mitosisTimeRangePercent = c(60,90);
gap = 4
ymax=5

#phase options / initialisation
palette = brewer.pal(7, "Set1")
colorcode = c(palette[1], palette[5], palette[3], palette[2], palette[4], "darkblue", "black")
includeDiffCells = FALSE;
phases = c("G1","S","G2","M","MeiS","Diff","Total");

indices = {}
indices[[1]]=c(2,6)
indices[[2]]=3
indices[[3]]=4
indices[[4]]=5
indices[[5]]=7
indices[[6]]=8
if(includeDiffCells){
	indices[[7]] = c(2,3,4,5,6,7,8);
}else{
	indices[[7]] = c(2,3,4,5,6,7);
}



par(mfrow=c(length(dirs),3))
for(d in dirs){

	testdata = read.table(paste(resultsDirectory,d,repNums[1],"/results_from_time_0/PhaseData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
	#testdata = read.table(paste(resultsDirectory,d,"/results_from_time_0/PhaseData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
    
    times = testdata[,1]
	tRange1 = times[floor((mitosisTimeRangePercent[1]/100)*length(times))]
	tRange2 = times[floor((mitosisTimeRangePercent[2]/100)*length(times))]
	tRange = seq(tRange1,tRange2,by=gap);
	total = data.frame()
	mitoses = data.frame()

    startingList = {}
    startingList[[1]] = rep(0,length(times))
    startingList[[2]] = rep(0,length(times))
    
    phaseTotals = {}
    for(i in seq(1,7,by=1)){
        phaseTotals[[i]] = startingList
    }
  

    #Loop over repeats ------------------

    for(rep in repNums){
    
        mitosisdata<-read.table(paste(resultsDirectory,d,rep,"/results_from_time_0/MitosisData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
        phasedata<-read.table(paste(resultsDirectory,d,rep,"/results_from_time_0/PhaseData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
 
        #mitosisdata<-read.table(paste(resultsDirectory,d,"/results_from_time_0/MitosisData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
        #phasedata<-read.table(paste(resultsDirectory,d,"/results_from_time_0/PhaseData.txt",sep=""),as.is=TRUE,header=FALSE,sep="\t");
 
        for(p in seq(1,7,by=1)){
 
        	sum = rep(0,length(times));
        	for(index in indices[[p]]){
        		sum = sum + phasedata[,index];
        	}
 
        	phaseTotals[[p]][[1]] = phaseTotals[[p]][[1]] + sum;
        	phaseTotals[[p]][[2]] = phaseTotals[[p]][[2]] + sum * sum;
 
        }

        
        for(t in seq(1,length(tRange)-1,by=1)){

        	selectedCells = mitosisdata[(mitosisdata[,1]>tRange[t] & mitosisdata[,1]<tRange[t+1]),]
        	selectedCells = selectedCells[2:length(selectedCells)]
        	nRows = (dim(selectedCells)[2])/2
        
        	if(rep==repNums[1]){
        		total = rbind( total, colSums(selectedCells[ 1:nRows ]) )
        		mitoses = rbind( mitoses, colSums(selectedCells[ (nRows+1):(2*nRows) ]) )

        	}else{
        		total[t,] = total[t,] + colSums(selectedCells[ 1:nRows ])
        		mitoses[t,] = mitoses[t,] + colSums(selectedCells[ (nRows+1):(2*nRows) ])
        	}

        }

    }

	plot(0,0, col="white", xlab="Cell row", ylab="Mitotic Index (%)", ylim=c(0,ymax), xlim=c(0,nRows), main="Mitotic profile")
	N = length(repNums)

	for(t in seq(1,length(tRange)-1,by=1)){
		
		c=rgb(t/length(tRange), 0, 1-t/length(tRange),1)
		MI = 100*( (mitoses[t,]/N) / (total[t,]/N) )
		points(seq(1,nRows,by=1), MI, col=c, type='l', lwd=2)
	}

	grad = rgb( seq(1,length(tRange)-1)/length(tRange), 0, 1-seq(1,length(tRange)-1)/length(tRange), 1);
	gradient.rect(40,ymax*0.5,50,ymax*0.9,col=grad,gradient="y")
	text(40,ymax*0.48,labels=paste("t=",tRange1))
	text(40,ymax*0.92,labels=paste("t=",tRange2))
	points(c(seq(1,10,by=1),15,20),100*c(0,4/367,4/545,8/592,7/644,12/660,12/667,11/645,7/660,8/643,4/643,0),type='l',lwd=2,col="black")


    #Plotting and processing ------------

    plot(0,0, main="Cell cycle phase percentages over time", xlab="Time (hph)", ylab="Percentage cells", ylim=c(0,100), xlim=c(min(times),max(times)), col="white")
    
    for(p in seq(1,7,by=1)){
    
    	phaseTotals[[p]][[1]] = phaseTotals[[p]][[1]]/length(repNums);
    	phaseTotals[[p]][[2]] = sqrt( phaseTotals[[p]][[2]]/length(repNums) - phaseTotals[[p]][[1]]^2);
    
    }

    for(p in seq(1,5,by=1)){
    
    	percentMean = 100*(phaseTotals[[p]][[1]]/phaseTotals[[7]][[1]]);
    	plusSD = percentMean + 100*(phaseTotals[[p]][[2]]/phaseTotals[[7]][[1]])
    	minusSD = percentMean - 100*(phaseTotals[[p]][[2]]/phaseTotals[[7]][[1]])
    	points(times, percentMean, type='l', lwd=2, col=colorcode[p]);
    	polygon(c(times,rev(times)), c(plusSD,rev(minusSD)), density = 40, angle=90, border=NA, col=colorcode[p], lwd=1)
    
    }
    
    if(includeDiffCells){
  
    	percentMean = 100*(phaseTotals[[6]][[1]]/phaseTotals[[7]][[1]]);
    	plusSD = percentMean + 100*(phaseTotals[[6]][[2]]/phaseTotals[[7]][[1]])
    	minusSD = percentMean - 100*(phaseTotals[[6]][[2]]/phaseTotals[[7]][[1]])
    	points(times, percentMean, type='l', lwd=2, col=colorcode[6]);
    	polygon(c(times,rev(times)), c(plusSD,rev(minusSD)), density = 40, angle=90, border=NA, col=colorcode[6], lwd=1)
    	legend("topright", legend = phases[1:6], fill = colorcode[1:6])

    }else{
    	legend("topright", legend = phases[1:5], fill = colorcode[1:5])
    }

    plot(times,phaseTotals[[7]][[1]],type='l',lwd=2,xlab="Time (hph)", ylab="Cell count", main = "Total proliferative cells");
    polygon(c(times,rev(times)), c(phaseTotals[[7]][[1]]+phaseTotals[[7]][[2]],rev(phaseTotals[[7]][[1]]-phaseTotals[[7]][[2]])), density = 40, angle=90, border=NA, col="black", lwd=1)
}


