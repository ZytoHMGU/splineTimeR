networkProperties <- function(igr) {
   
   if(!xor(is(igr, "igraph"), all(lapply(igr, class) == "igraph"))) stop("igr must be an igraph object or list of igraph objects")
   
   if(is(igr,"list") & is.null(names(igr))) {
      names(igr) = paste("igraph_", 1:length(igr), sep="")
   }

   temp <- c()
   temp2 <- c()
   
   if(is(igr,"list")) {
      for(i in 1:length(igr)) {
         igr.temp = assign(paste("igr.", i, sep=""), igr[[i]])
         temp <- c(temp, length(V(igr.temp)$name))
         temp2 <- c(temp2, length(E(igr.temp)))
      }
   } else {
         igr.1 = igr
         temp <- c(temp, length(V(igr.1)$name))
         temp2 <- c(temp2, length(E(igr.1)))
   }
      
   ## load Functional Interaction data for BioGRID and Reactome
   data(FIs)

   ## Reactome
   reactome <- FIs$FIs_Reactome
   igr_reactome <- graph.data.frame(reactome, directed=FALSE, vertices=NULL)
   geneNumber_reactome <- length(unique(c(as.character(reactome$V1), as.character(reactome$V2))))
   degree_reactome <- degree(igr_reactome, v=V(igr_reactome), mode="all", loops=FALSE, normalized=FALSE)
   freq_of_deg_reactome <- table(degree_reactome)
   
   ## bioGrid
   bioGrid <- FIs$FIs_BioGRID
   igr_bioGrid <- graph.data.frame(bioGrid, directed=FALSE, vertices=NULL)
   geneNumber_bioGrid <- length(unique(c(as.character(bioGrid$V1), as.character(bioGrid$V2))))
   degree_bioGrid <- degree(igr_bioGrid, v=V(igr_bioGrid), mode="all", loops=FALSE, normalized=FALSE)
   freq_of_deg_bioGrid <- table(degree_bioGrid)
   
   temp <- c(temp, length(colnames(igr_reactome[,])), length(colnames(igr_bioGrid[,])))
   temp2 <- c(temp2,length(E(igr_reactome)), length(E(igr_bioGrid)))
   
   ## user data
   if(is(igr,"list")) {
      for(i in 1:length(igr)) {
         igr.temp <- get(paste("igr.",i,sep=""))
         degree_temp <- assign(paste("degree.",i,sep=""), degree(igr.temp, v=V(igr.temp), mode="all", loops=FALSE, normalized=FALSE))
         assign(paste("freq_of_deg.",i,sep=""), table(degree_temp))
      }
   } else {
      degree.1 <- degree(igr.1, v=V(igr.1), mode="all", loops=FALSE, normalized=FALSE)
      freq_of_deg.1 <- table(degree.1)
   }

   #### degree distribution ####
   pdf(paste("degree_distribution.pdf",sep=""), width=12, height=6)
   par(mfrow = c(1,2))
   
   if(is(igr,"list")) {
      for(i in 1:length(igr)) {
         plot(ecdf(get(paste("degree.",i,sep=""))), main=paste("degree distribution - ecdf\n", names(igr)[i], sep=""))
         hist(get(paste("degree.",i,sep="")), main=paste("degree distribution - histogram\n", names(igr)[i], sep=""), breaks=50, xlab="degree")
      }
   } else {
      plot(ecdf(degree.1), main=paste("degree distribution - ecdf\n igraph",sep=""))
      hist(degree.1, main=paste("degree distribution - histogram\n igraph", sep=""), breaks=50, xlab="degree")
   }
   
   ## Reactome	
   plot(ecdf(degree_reactome), main="degree distribution - ecdf\nReactom")
   hist(degree_reactome, main="degree distribution - histogram \nReactom", breaks=50, xlab="degree")
   
   dev.off()
   
   #### scale-freeness ####
   temp3 = c()
   pdf(paste("scale_free_properties.pdf", sep=""), width=6.5, height=6.5)
   
   ## user network(s)
   if(is(igr,"list")) {
      for(i in 1:length(igr)) {
         plot(as.numeric(names(get(paste("freq_of_deg.",i,sep="")))), get(paste("freq_of_deg.",i,sep=""))/sum(get(paste("freq_of_deg.",i,sep=""))), xlab="k", ylab="P(k)", main=paste("degree distribution\n", names(igr)[i], sep=""), log="xy", yaxt="n", pch=20)
         lm_fit = lm(log10(get(paste("freq_of_deg.",i,sep=""))/sum(get(paste("freq_of_deg.",i,sep="")))) ~ log10(as.numeric(names(get(paste("freq_of_deg.",i,sep=""))))))
         abline(lm_fit)
         legend("topright", bty="n", legend=paste("R2 = ", format(summary(lm_fit)$adj.r.squared, digits=4),"\n y = ", format(lm_fit$coefficients[2], digits=4),"x + ", format(lm_fit$coefficients[1], digits=4), sep=""))
         axis(2, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
         
         temp3 <- c(temp3, -round(lm_fit$coefficients[2], digits=4))
      }
   } else {
      plot(as.numeric(names(freq_of_deg.1)), freq_of_deg.1/sum(freq_of_deg.1), xlab="k", ylab="P(k)", main=paste("degree distribution\n igraph", sep=""), log="xy", yaxt="n", pch=20)
      lm_fit = lm(log10(freq_of_deg.1/sum(freq_of_deg.1)) ~ log10(as.numeric(names(freq_of_deg.1))))
      abline(lm_fit)
      legend("topright", bty="n", legend=paste("R2 = ", format(summary(lm_fit)$adj.r.squared, digits=4),"\n y = ", format(lm_fit$coefficients[2], digits=4),"x + ", format(lm_fit$coefficients[1], digits=4), sep=""))
      axis(2, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
      
      temp3 <- c(temp3, -round(lm_fit$coefficients[2], digits=4))
   }
   
   ## Reactome
   plot(as.numeric(names(freq_of_deg_reactome)), freq_of_deg_reactome/sum(freq_of_deg_reactome), xlab="k", ylab="P(k)", main="degree distribution\nReactom", log="xy", yaxt="n", pch=20)
   lm_fit = lm((log10(freq_of_deg_reactome/sum(freq_of_deg_reactome))) ~ log10(as.numeric(names(freq_of_deg_reactome))))
   abline(lm_fit)
   legend("topright", bty="n", legend=paste("R2 = ", format(summary(lm_fit)$adj.r.squared, digits=4),"\n y = ", format(lm_fit$coefficients[2], digits=4),"x + ", format(lm_fit$coefficients[1], digits=4), sep=""))
   axis(2, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
   
   temp3 <- c(temp3, -round(lm_fit$coefficients[2], digits=4))
   
   ## BioGRID
   plot(as.numeric(names(freq_of_deg_bioGrid)), freq_of_deg_bioGrid/sum(freq_of_deg_bioGrid), xlab="k", ylab="P(k)", main="degree distribution\nBioGRID", log="xy", yaxt="n", pch=20)
   lm_fit = lm((log10(freq_of_deg_bioGrid/sum(freq_of_deg_bioGrid))) ~ log10(as.numeric(names(freq_of_deg_bioGrid))))
   abline(lm_fit)
   legend("topright", bty="n", legend=paste("R2 = ", format(summary(lm_fit)$adj.r.squared, digits=4),"\n y = ", format(lm_fit$coefficients[2], digits=4),"x + ", format(lm_fit$coefficients[1], digits=4), sep=""))
   axis(2, at = c(1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1))
   
   temp3 <- c(temp3, -round(lm_fit$coefficients[2], digits=4))
   dev.off()
 
   if(is(igr,"list")) {  
      info.temp <- cbind(temp,temp2,temp3)
      rownames(info.temp) <- c(names(igr),"Reactome","BioGRID")
   } else {
      info.temp <- cbind(temp,temp2,temp3)
      rownames(info.temp) <- c("igraph","Reactome","BioGRID")
   }

   colnames(info.temp) <- c("nodes","edges","degree_exponent")
   return(info.temp)
}