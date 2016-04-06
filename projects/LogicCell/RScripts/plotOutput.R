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
par( mfrow=c(ceiling(sqrt(N)),ceiling(sqrt(N))), mar=c(1,1,1,1), mgp=c(2,0.3,0) )

slidingWindowAv = FALSE
windowRad = 10 


#Processing and plotting loop
for( item in seq(1,N,by=1) ){

	outdata<-read.table(paste(resultsDirectory, dirs[item], "/results_from_time_0/OutputData.txt",sep=""),
		     as.is=TRUE,header=FALSE,sep="\t");

	spl = strsplit(dirs[item],"_")
	z = spl[[1]][2] 
	s = spl[[1]][3] 
	r = spl[[1]][4] 
	header = paste("Resist: ", r, " Zone: ", z, " Spring: ", s)


	if(!slidingWindowAv){

		plot(outdata[,1], outdata[,2], 
			 type='l', xlab="Time (hph)", ylab="Output (cells per hour)",lwd=2,main=header);

		points(c(1,max(outdata[,1])),c(15.5,15.5),type='l',col='red',lty=3);	

	}else{

		data = rep(0, length(outdata[,1]))

		for( t in seq(windowRad+1,length(data)-windowRad-1,by=1) ){

			data[t] = sum(outdata[seq(t-windowRad,t+windowRad,by=1),2])/(2*windowRad);
		}

		plot(outdata[seq(windowRad+1,length(data)-windowRad-1,by=1),1], data[seq(windowRad+1,length(data)-windowRad-1,by=1)],
		     type='l', xlab="Time (hph)", ylab="Output (cells per hour)",lwd=2, main=header)
	
		points(c(windowRad+1,length(data)-windowRad-1),c(15.5,15.5),type='l',col='red',lty=3);	
	}
}