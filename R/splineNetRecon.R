splineNetRecon <- function(eSetObject, treatmentType, probesForNR="all", cutoff.ggm=0.8, method="dynamic", saveEdges=FALSE) {

	if(!is(eSetObject, "ExpressionSet")) stop("eSetObject must be of class ExpressionSet")
   if(!(("SampleName" %in% names(pData(eSetObject))) & ("Time" %in% names(pData(eSetObject))) & ("Treatment" %in% names(pData(eSetObject))))) stop("eSetObject has to include SampleName, Time and Treatment columns in phenotypic data")
	if(!("Replicate" %in% colnames(pData(eSetObject)))) {
	   pData(eSetObject)$Replicate = rep("A",nrow(pData(eSetObject)))
		cat("\n----------------------------------------------------------------\n")
		cat("No replicates identified","\n")
	}

	if(all(probesForNR == "all")) {
		probesForNR = rownames(exprs(eSetObject))
		cat("\n----------------------------------------------------------------\n")
		cat("All rows from eSetObject will be taken for network reconstruction\n")
		
	} else {
		if(!is(probesForNR, "character")) stop("Define valid rownames of exprs(eSetObject) to plot")
	}
   
	if(!all(probesForNR %in% rownames(exprs(eSetObject)))) stop("Some of provided rownames don't exist in eSetObject")
   if(!is(cutoff.ggm,"numeric")) stop("cutoff.ggm must be numeric")
	if(!(treatmentType %in% unique(pData(eSetObject)$Treatment))) stop("Choose valid treatment")
	if(!(method %in% c("dynamic","static"))) stop("Define method for network reconstruction: \"dynamic\" or \"static\"")
   if(!is(saveEdges,"logical")) stop("saveEdges must be logical")

### prepare data ###
	# introduce new name for ordering samples according to time and replicates
	colnames(exprs(eSetObject)) <- paste("T", pData(eSetObject)$Time, "_", pData(eSetObject)$Treatment, "_", pData(eSetObject)$Replicate, sep="")
	exprs.NR = exprs(eSetObject)[probesForNR, which(pData(eSetObject)$Treatment == treatmentType)]
	exprs.NR = exprs.NR[,mixedsort(colnames(exprs.NR))]

	# introduce new help variable for identifying numbers of replicates 
	help_repeats <- paste("T", pData(eSetObject)$Time, "_", pData(eSetObject)$Treatment, sep="")
	help_repeats = help_repeats[which(pData(eSetObject)$Treatment == treatmentType)]
	help_repeats = mixedsort(help_repeats)

	# check for replicates
	repeats = c()
	for(i in unique(help_repeats)) {
		repeats = c(repeats, length(which(help_repeats == i)))
	}

	exprs.NR = t(exprs.NR)
	node.labels.NR = probesForNR
	
### create longitudinal object ###
	time = unique(pData(eSetObject)$Time[which(pData(eSetObject)$Treatment == treatmentType)])
	time = time[order(time)]

	exprs.NR = as.longitudinal(exprs.NR, repeats=repeats, time = time)
	
	cat("\n----------------------------------------------------------------\n")
	cat("Longitudinal object")
	cat("\n----------------------------------------------------------------\n")
	print(get.time.repeats(exprs.NR))
	
### estimate partial correlation matrix ###
	pcor.NR <- ggm.estimate.pcor(exprs.NR, method=method)

### test edges ###
	edges.NR <- network.test.edges(pcor.NR, fdr=TRUE, direct=FALSE, plot=FALSE)
	if(saveEdges) {
	   save(edges.NR, node.labels.NR, file=paste("edge_list_", method, ".Rdata", sep=""))
	}

### create igraph object ###	
	if(length(cutoff.ggm)!=1) {
	   igr <- list()
	   j <- 1
	   for(i in cutoff.ggm) {
	      cat(paste("\n-------------- ","igraph_",i," --------------",sep=""))
	      net.NR <- extract.network(edges.NR, cutoff.ggm=i)
	      igr.help <- network.make.graph(edge.list=net.NR, node.labels=node.labels.NR, drop.singles=TRUE)
	      igr.help <- igraph.from.graphNEL(graphNEL=igr.help, name=TRUE, weight=FALSE, unlist.attrs=TRUE)
	      cat("Number of nodes: ", length(colnames(igr.help[,])),"\n")
	      igr[[j]] <- igr.help
	      names(igr)[j] <- paste("igraph_",i,sep="")
	      j=j+1
	   }
	} else {
	   cat(paste("\n-------------- igraph_",cutoff.ggm," --------------",sep=""))
	   net.NR <- extract.network(edges.NR, cutoff.ggm=cutoff.ggm)
	   igr <- network.make.graph(edge.list=net.NR, node.labels=node.labels.NR, drop.singles=TRUE)
	   igr <- igraph.from.graphNEL(graphNEL=igr, name=TRUE, weight=FALSE, unlist.attrs=TRUE)
	   cat("Number of nodes: ", length(colnames(igr[,])),"\n")
	}
	return(igr)
}