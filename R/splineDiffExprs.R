splineDiffExprs <- function(eSetObject, df, cutoff.adj.pVal=1, reference, intercept=TRUE) { 
	
	if(!is(eSetObject, "ExpressionSet")) stop("eSetObject must be of class ExpressionSet")
   if(!(("SampleName" %in% names(pData(eSetObject))) & ("Time" %in% names(pData(eSetObject))) & ("Treatment" %in% names(pData(eSetObject))))) stop("eSetObject has to include SampleName, Time and Treatment columns in phenotypic data")
	if(!(is(cutoff.adj.pVal, "numeric") & (cutoff.adj.pVal >= 0) & (cutoff.adj.pVal <= 1))) stop("cutoff.adj.pVal must be numeric between 0 and 1")
	if(!(is(df, "numeric") & (df%%1 == 0) & (df > 0))) stop("df must be integer > 0")
   if(!(reference %in% levels(factor(pData(eSetObject)$Treatment)))) stop("define valid reference group")
	if(!is(intercept, "logical")) stop("intercept must be boolean data type")
	
	b_ <- ns(pData(eSetObject)$Time, df=df)
	d_ <- factor(pData(eSetObject)$Treatment, levels=c(reference, setdiff(levels(factor(pData(eSetObject)$Treatment)), reference)))
	design <- model.matrix(~d_*b_)
	fit <- lmFit(eSetObject, design)
	fit_eBayes <- eBayes(fit)	
	
	if(ncol(fData(eSetObject))==0) {
	   fData(eSetObject) <- data.frame(rownames(exprs(eSetObject)))
	   colnames(fData(eSetObject)) <- "row_IDs" 
	}

	if(intercept) {
	   fit_coeff_ref <- fit$coefficient[,c(1,3:(df+2))]
	   colnames(fit_coeff_ref)[1] <- "b_0"
	   fData(eSetObject) <- cbind(fData(eSetObject), fit_coeff_ref)
	   topTable_output <- topTable(fit_eBayes, coef=c(2,(df+3):(2*df+2)), sort.by="none", number=Inf, genelist=fData(eSetObject))
		colnames(topTable_output)[(ncol(topTable_output)-4-df):(ncol(topTable_output)-4)] <- paste("d_",0:df,sep="")
	} else {
	   fit_coeff_ref <- fit$coefficient[,c(1,3:(df+2),2)]
	   colnames(fit_coeff_ref)[1] <- "b_0"
	   fData(eSetObject) <- cbind(fData(eSetObject), fit_coeff_ref)
	   topTable_output <- topTable(fit_eBayes, coef=c((df+3):(2*df+2)), sort.by="none", number=Inf, genelist=fData(eSetObject))
	   colnames(topTable_output)[(ncol(topTable_output)-4-df):(ncol(topTable_output)-4)] <- paste("d_",0:df,sep="")
	}
	
	diffExprs <- topTable_output[topTable_output$adj.P.Val <= cutoff.adj.pVal,]
	diffExprs <- diffExprs[order(diffExprs$adj.P.Val),]
	cat("-------------------------------------------------", "\n")
 	cat(paste("Differential analysis done for df = ", df," and adj.P.Val <= ", cutoff.adj.pVal, sep=""), "\n")
 	cat("Number of differentially expressed genes: ", nrow(diffExprs),"\n")

	return(diffExprs)
}