powers = seq(14,6,-1)
population = "node"
cycle = "stoch"
method = "RK4"
adaptive = "false"
reps = seq(0,4)

#-----------------------------------------------------------------------

dir = "/home/kathryn/ChasteFork/testoutput/"

for(p in powers){

	timings = c()

	for(r in reps){
	
		file = paste(population,"_",cycle,"_",method,"_",adaptive,"_p",p,"_",r,sep="")
		path = paste(dir,file,"/results_from_time_0", sep="")
		allVTUFiles = list.files(path, pattern="*.vtu")
		fileNumbers = as.numeric(gsub("[^0-9]", "",allVTUFiles))
		maxval = max(fileNumbers)

		pathmin = paste(path,"/results_0.vtu", sep="")
		pathmax = paste(path,"/results_",maxval,".vtu", sep="")
		
		minTime = file.info(pathmin)$ctime
		maxTime = file.info(pathmax)$ctime

		diffMins = difftime(maxTime, minTime, units="mins")
		timings = c(timings,diffMins)
	}
	
	mymax = max(timings)
	mymean = mean(timings)
	#TODO: outlier detection

	print(paste("Power: ", p))
	print(paste("Average run time: ", mymean))
	print(paste("Max run time: ", mymax))
}

	