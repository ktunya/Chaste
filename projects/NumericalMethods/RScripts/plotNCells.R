power = 14
population = "node"
cycle = "stoch"
method = "FE"
reps = seq(0,4)
average = TRUE

#-----------------------------------------------------------------------

dir = "/home/kathryn/ChasteFork/testoutput/"
powers = seq(14,6,-1)
mycolors = c(rgb(0.36,0.56,0.99,1.0),rgb(0.26,0.87,1.0,1.0),rgb(0.0,1.0,0.46,1.0),
	         rgb(0.0,0.64,0.0,1.0),rgb(0.77,1.0,0.0,1.0),rgb(0.96,0.78,0.0,1.0),
	         rgb(0.98,0.29,0.016,1.0),rgb(0.98,0.0,0.28,1.0),rgb(0.79,0.0,0.59,1.0));
color = rev(mycolors)[powers==power]
mylabel = c(paste("2^",power,sep=""));

if(average == TRUE){

	for(r in reps){

		file = paste(population,"_",cycle,"_",method,"_false_p",power,"_",r,sep="")
		path = paste(dir,file,"/results_from_time_0/PositionData.txt", sep="")
		data = read.table(path, as.is=TRUE, header=FALSE, sep="\t")

		if(r==0){
			times = unique(data[,1])
			nCells = rep(0, length(times))
		}

		index = 1
		for(t in times){
			cellmax = max(data[data[,1]==t,2]) 
			nCells[index] = nCells[index] + cellmax
			index = index+1
		}
	}

	nCells = nCells / length(reps)

}else{
	file = paste(population,"_",cycle,"_",method,"_false_p",power,"_0", sep="")
	path = paste(dir,file,"/results_from_time_0/PositionData.txt", sep="")
	data = read.table(path, as.is=TRUE, header=FALSE, sep="\t")

	times = unique(data[,1])
	nCells = rep(0, length(times))
	index = 1
	for(t in times){
		cellmax = max(data[data[,1]==t,2]) 
		nCells[index] = cellmax
		index = index+1
	}
}	

plot(times, nCells, lwd=2, type='l', col=color,
     xlab="Time (hours)", ylab="Number of cells")
	
legend("topleft", fill=color, legend=mylabel)