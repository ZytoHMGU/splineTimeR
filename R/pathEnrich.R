pathEnrich <- function(geneList, geneSets, universe=NULL) {
   
   if(!is(geneList, "character")) stop("geneList must be a vector of gene names")
   if(!is(geneSets, "GeneSetCollection")) stop("geneSets must be of class GeneSetCollection")
   if(!(is(universe, "numeric") | is.null(universe))) stop("universe must be numeric")
   
   info <- c()
   for(i in c(1:length(geneSets))) info[i] <- description(geneSets[[i]])
   geneSets <- geneIds(geneSets)
   
   
   genesInGeneSets <- unique(unlist(geneSets))
   geneList <- as.character(geneList[geneList %in% genesInGeneSets])
   
   if(is.null(universe)) universe <- length(genesInGeneSets)
   
   cat("--------------------------------------------------------","\n")
   
   genes_in_pathway <- lapply(geneSets, function(x) length(x))
   matches <- lapply(geneSets, function(x) length(which(!is.na(match(x,geneList)))))
   index <- lapply(geneSets, function(x) which(!is.na(match(x,geneList))))
   
   overlap = c() 
   for(i in c(1:length(index))) {
      if(length(index[[i]]) != 0) {
         overlap[i] <- paste(geneSets[[i]][index[[i]]], collapse=",")
      } else {
         overlap[i] <- ""
      }
   }
   
   percent <- round((unlist(matches)/unlist(genes_in_pathway))*100, 2)
   input <- cbind(unlist(genes_in_pathway), unlist(matches))
   pValue <- do.call("phyper", list(q=input[,2]-1, m=length(geneList), n=universe-length(geneList), k=input[,1], lower.tail=FALSE, log.p=FALSE), quote=FALSE)
   
   # pValue correction
   adj.pValue <- p.adjust(pValue, method = "BH")
   
   # construct result table
   enrichPath <- cbind(rownames(input), info, input, percent, pValue, adj.pValue, overlap)
   enrichPath <- data.frame(enrichPath, stringsAsFactors = FALSE)
   colnames(enrichPath) <- c("pathway", "description", "genes_in_pathway", "matches", "%_match", "pValue", "adj.pValue", "overlap")
   
   cat("Pathway enrichment done!","\n")
   cat("--------------------------------------------------------","\n")
   enrichPath <- enrichPath[with(enrichPath, order(adj.pValue)),]
   
   write.table(enrichPath, file="pathway_enrichment_result.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE) 
   
   return(enrichPath)
}