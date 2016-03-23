\name{pathEnrich}
\alias{pathEnrich}
\title{Pathway enrichment analysis}

\description{
Function performs a pathway enrichment analysis of a definied set of genes.
}
\usage{pathEnrich(geneList, geneSets, universe=NULL)}

\arguments{
  \item{geneList}{vector of gene names to be used for pathway enrichment}
  \item{geneSets}{\code{"GeneSetColletion"} object with functional pathways gene sets}
  \item{universe}{number of genes that were probed in the initial experiment}
} 

\details{
\code{geneSets} is a \code{"GeneSetColletion"} object containing gene sets from various databases. Different sources for gene sets data are allowed and have to be provided in Gene Matrix Transposed file format (*.gmt), where each gene set is described by a pathway name, a description, and the genes in the gene set. Two examples are shown to demonstrate how to define \code{geneSets} object. See examples.

The variable \code{universe} represents a total number of genes that were probed in the initial experiment, e.g. the number of all genes on a microarray. If \code{universe} is not definied, \code{universe} is equal to the number of all genes that can be mapped to any pathways in chosen database.
}

\value{
A data.frame with following columns:
\item{pathway}{names of enriched pathways}
\item{description}{gene set description (e.g. a link to the named gene set in MSigDB)}
\item{genes_in_pathway}{total number of known genes in the pathway}
\item{\%_match}{number of matched genes refered to the total number of known genes in the pathway given in \%}
\item{pValue}{p-value}
\item{adj.pValue}{Benjamini-Hochberg adjucted p-value}
\item{overlap}{genes from input genes list that overlap with all known genes in the pathway}

Additionally an .txt file containing all above information is created.
}

\author{Agata Michna}

\references{
Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., Paulovich, A., Pomeroy, S. L., Golub, T. R., Lander, E. S. and Mesirov, J. P. (2005). Gene set enrichment analysis: A knowledge-based approach for interpreting genome-wide expression profiles.
\emph{PNAS} 102(43), 15545-15550.

\url{http://www.broadinstitute.org/gsea/msigdb/collections.jsp}

\url{http://www.reactome.org/pages/download-data/}
}

\examples{
\dontrun{
   ## Example 1 - using gene sets from the Molecular Signatures Database (MSigDB collections)
   ## Download .gmt file 'c2.all.v5.0.symbols.gmt' (all curated gene sets, gene symbols)
   ## from the Broad, http://www.broad.mit.edu/gsea/downloads.jsp#msigdb, then
   geneSets <- getGmt("/path/to/c2.all.v5.0.symbols.gmt")
   ## load "eSetObject" containing simulated time-course data
   data(TCsimData)
   ## check for differentially expressed genes
   diffExprs <- splineDiffExprs(eSetObject = TCsimData, df = 3, cutoff.adj.pVal = 0.01, reference = "T1")
   ## use differentially expressed genes for pathway enrichment analysis
   enrichPath <- pathEnrich(geneList = rownames(diffExprs), geneSets = geneSets, universe = 6536)}

\dontrun{
   ## Example 2 - using gene sets from the Reactome Pathway Database
   ## Download and unzip .gmt.zip file 'ReactomePathways.gmt.zip'
   ## ("Reactome Pathways Gene Set" under "Specialized data formats") from the Reactome website
   ## http://www.reactome.org/pages/download-data/, then
   geneSets <- getGmt("/path/to/ReactomePathways.gmt")
   data(TCsimData)
   diffExprs <- splineDiffExprs(eSetObject = TCsimData, df = 3, cutoff.adj.pVal = 0.01, reference = "T1")
   enrichPath <- pathEnrich(geneList = rownames(diffExprs), geneSets = geneSets, universe = 6536)}
   
## Small example with gene sets consist of KEGG pathways only
geneSets <- getGmt(system.file("extdata", "c2.cp.kegg.v5.0.symbols.gmt", package="splineTCDiffExpr"))
data(TCsimData)
diffExprs <- splineDiffExprs(eSetObject = TCsimData, df = 3, cutoff.adj.pVal = 0.01, reference = "T1")
enrichPath <- pathEnrich(geneList = rownames(diffExprs), geneSets = geneSets, universe = 6536)
}

\keyword{gene set enrichment analysis}
\keyword{pathway enrichment analysis}
