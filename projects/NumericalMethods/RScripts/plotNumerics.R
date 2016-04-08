powers = seq(6,8,1)
population = "mesh"
cycle = "stoch"
method = "AM2"
adaptive = "false"
nReps = 1
norm = "LInf"
save = FALSE
#-------------------------------------------------------------------------------

positionColumns = c(3,4)
if(population=="node" || population=="particle"){
	positionColumns = c(3,4,5)
}

name = paste(population,"_",cycle,"_",method,"_",adaptive,".eps")
reps = seq(0,nReps-1,1)
dir = "/home/kathryn/ChasteFork/testoutput/"

mycolors = rev(c(rgb(0.36,0.56,0.99,1.0),rgb(0.26,0.87,1.0,1.0),rgb(0.0,1.0,0.46,1.0),
	         rgb(0.0,0.64,0.0,1.0),rgb(0.77,1.0,0.0,1.0),rgb(0.96,0.78,0.0,1.0),
	         rgb(0.98,0.29,0.016,1.0),rgb(0.98,0.0,0.28,1.0),rgb(0.79,0.0,0.59,1.0)));
mylabels=c()
for(p in powers){
	mylabels = c(mylabels,paste("2^",p,sep=""));
}
mylabels = mylabels

if(save){
	postscript(name, horiz=TRUE, bg="white", paper="special", width=4.5, height=4.5, pointsize=12)
}


for(r in reps){
	
	baselineFile = paste(population, "_", cycle, "_RK4_false_p12_", r, sep="")
	print(baselineFile)
	baselinePath = paste(dir,baselineFile,"/results_from_time_0/PositionData.txt", sep="")
	baseline = read.table(baselinePath, as.is=TRUE, header=FALSE, sep="\t")

	if(r==0){
		baselineTimes = unique(baseline[,1])
		errorMat = matrix(0, length(powers), length(baselineTimes))
	}

	pindex = 1
	for(pow in powers){
	
		file = paste(population, "_", cycle, "_", method, "_", adaptive, "_p", pow, "_", r, sep="")
		fullPath = paste(dir, file, "/results_from_time_0/PositionData.txt", sep="")
	
		print(file)
	
		data = read.table(fullPath, as.is=TRUE, header=FALSE, sep="\t")
		times = unique(data[,1])
	
		tindex = 1
		for(t in times){
			cells = intersect(unique(data[data[,1]==t,2]), unique(baseline[baseline[,1]==t,2]))
			positionsOld = baseline[(baseline[,1]==t & baseline[,2]%in%cells), positionColumns]
			positionsNew = data[(data[,1]==t & data[,2]%in%cells), positionColumns]
			diff = positionsOld - positionsNew
			sqDisp = apply(diff*diff,1,sum)
			if(norm=="LInf"){
				maxDiff = max(sqrt(sqDisp))
				errorMat[pindex,tindex] = errorMat[pindex,tindex] + maxDiff
			}else if(norm=="L2"){
				meanDiff = mean(sqrt(sqDisp))
				errorMat[pindex,tindex] = errorMat[pindex,tindex] + meanDiff
			}else{
				print("Unrecognised norm")
			}
			tindex = tindex+1
		}
		
		pindex = pindex+1	
	}
}
errorMat = errorMat/length(reps)


for(p in seq(1,length(powers))){

	error = errorMat[p,]

	if(powers[p] == powers[1]){
		plot(baselineTimes, error, lwd=2, type='l', col=mycolors[p], ylim=c(0,3),
		     xlab="Time (hours)", ylab=paste(norm," norm of cell position error")) 
	}else{
		points(baselineTimes, error, lwd=2, type='l', col=mycolors[p])
	}
}	
legend("topleft", fill=mycolors, legend=mylabels)

if(save){
	dev.off()
}