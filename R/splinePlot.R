splinePlot <- function(eSetObject, df, reference, toPlot="all") {

	if(!is(eSetObject, "ExpressionSet")) stop("eSetObject must be of class ExpressionSet")
	if(!(("SampleName" %in% names(pData(eSetObject))) & ("Time" %in% names(pData(eSetObject))) & ("Treatment" %in% names(pData(eSetObject))))) stop("eSetObject has to include SampleName, Time and Treatment columns in phenotypic data")
	if(!(is(df, "numeric") & (df%%1 == 0) & (df > 0))) stop("df must be integer > 0")
   if(!(reference %in% levels(factor(pData(eSetObject)$Treatment)))) stop("define valid reference group")
   
	if(all(toPlot == "all")) {
		toPlot = rownames(exprs(eSetObject))
	} else {
		if(!is(toPlot, "character")) stop("define row names of exprs(eSetObject) to plot")
	}
	if(!all(toPlot %in% rownames(exprs(eSetObject)))) stop("some of provided names for plotting are not included in eSetObject")

	b_ <- ns(pData(eSetObject)$Time, df=df)
	d_ <- factor(pData(eSetObject)$Treatment, levels=c(reference, setdiff(levels(factor(pData(eSetObject)$Treatment)), reference)))
	design <- model.matrix(~d_*b_)
	fit <- lmFit(eSetObject, design)

	exprs.data <- exprs(eSetObject)
	factorTreatment <- levels(d_)
	
	timePoints_C <- unique(pData(eSetObject)$Time[pData(eSetObject)$Treatment == factorTreatment[1]])
	timePoints_T <- unique(pData(eSetObject)$Time[pData(eSetObject)$Treatment == factorTreatment[2]])

	regressionMatrix_C <- ns(timePoints_C, df=df)  
	regressionMatrix_T <- ns(timePoints_T, df=df)  

	newTime <- seq(min(c(timePoints_C,timePoints_T)), max(c(timePoints_C,timePoints_T)), length.out=101)

	regressionMatrixEval_C <- predict(regressionMatrix_C, newTime) 
	regressionMatrixEval_T <- predict(regressionMatrix_T, newTime)

	number = length(toPlot)
	legendComp = c(factorTreatment[1],factorTreatment[2])
	ylim = c(min(exprs.data[toPlot,])-0.25, max(exprs.data[toPlot,])+0.25)
	
	pdf(paste("plots_df",df,"_spline.pdf",sep=""), width=6.5, height=6.5)
	for(i in 1:number)
	{
		ix <- which(toPlot[i] == row.names(exprs.data))
		data_C <- exprs.data[ix,pData(eSetObject)$Treatment == factorTreatment[1]]
		data_T <- exprs.data[ix,pData(eSetObject)$Treatment == factorTreatment[2]]
		timePoints_C = pData(eSetObject)$Time[pData(eSetObject)$Treatment == factorTreatment[1]]
		timePoints_T = pData(eSetObject)$Time[pData(eSetObject)$Treatment == factorTreatment[2]]

		plot(timePoints_C, data_C, ylim=ylim, col=4, pch=20, main=paste(toPlot[i], sep="\n"), xlab="time", ylab="expression")
		points(timePoints_T, data_T, col=2, pch=20)
		legend("topright", lty=c(1,1), lwd=c(1.5,1.5), legendComp, col=c(4,2))

		coeffs <- fit$coefficient[ix,]
		newY_C <- coeffs[1]
		newY_T <- coeffs[1]+coeffs[2]

		for(i in c(3:(df*2+2-df))){
			newY_C <- newY_C + coeffs[i]*regressionMatrixEval_C[,(i-df+1)]
			newY_T <- newY_T + (coeffs[i]+coeffs[i+df])*regressionMatrixEval_T[,(i-df+1)]
		}
		lines(newTime, newY_C, col=4)
		lines(newTime, newY_T, col=2)
	}
	invisible(dev.off())
}