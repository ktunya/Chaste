powers = seq(14,6,-1)
population = "node"
cycle = "stoch"
method = "FE"
adaptive = "false"
nReps = 3

reps = seq(0,nReps-1,1)
dir = "/home/kathryn/ChasteFork/testoutput/"

mycolors = c(rgb(0.36,0.56,0.99,1.0),rgb(0.26,0.87,1.0,1.0),rgb(0.0,1.0,0.46,1.0),
	         rgb(0.0,0.64,0.0,1.0),rgb(0.77,1.0,0.0,1.0),rgb(0.96,0.78,0.0,1.0),
	         rgb(0.98,0.29,0.016,1.0),rgb(0.98,0.0,0.28,1.0),rgb(0.79,0.0,0.59,1.0));
mylabels=c()
for(p in powers){
	mylabels = c(mylabels,paste("2^",p,sep=""));
}

file = paste(population, "_", cycle, "_", method, "_", adaptive, "_p", powers[1], "_", nReps[1], sep="")
firstPath = paste(dir, file, "/results_from_time_0/PositionData.txt", sep="")
firstData = read.table(fullPath, as.is=TRUE, header=FALSE, sep="\t")
firstTimes = unique(firstData[,1])
errorMat = matrix(0, length(powers), length(firstTimes))


for(rep in reps){
	pindex = 1
	for(pow in powers){
	
		file = paste(population, "_", cycle, "_", method, "_", adaptive, "_p", pow, "_", rep, sep="")
		fullPath = paste(dir, file, "/results_from_time_0/PositionData.txt", sep="")
	
		print(file)
	
		data = read.table(fullPath, as.is=TRUE, header=FALSE, sep="\t")
		times = unique(data[,1])
	
		if(pow == powers[1]){
			baseline = data
		}else{
			tindex = 1
			for(t in times){
				cells = intersect(unique(data[data[,1]==t,2]), unique(baseline[baseline[,1]==t,2]))
				sqDiff = 0
				for(c in cells){
					positionOld = baseline[(baseline[,1]==t & baseline[,2]==c), c(3,4,5)]
					positionNew = data[(data[,1]==t & data[,2]==c), c(3,4,5)]
					diff = positionOld - positionNew
					sqDiff = sqDiff + sum(diff*diff)
				}
				errorMat[pindex,t] = errorMat[pindex,tindex] + sqrt(sqDiff)
				tindex = tindex+1
			}
		}
		pindex = pindex+1	
	}
}
errorMat = errorMat/length(reps)


for(p in seq(1,length(powers))){

	error = errorMat[p,]

	if(powers[p] == powers[1]){
		plot(firstTimes, error, lwd=2, type='l', col=mycolors[p], ylim=c(0,10),
		     xlab="Time (hours)", ylab="L2 norm error in cell position") 
	}else{
		points(firstTimes, error, lwd=2, type='l', col=mycolors[p])
	}
}	
legend("topleft", fill=rev(mycolors), legend=rev(mylabels))