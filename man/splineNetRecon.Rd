\name{splineNetRecon}
\alias{splineNetRecon}
\title{Network reconstruction based on partial correlation method with shrinkage approach}

\description{
\code{splineNetRecon} reconstructs gene association networks from time-course data. Based on given object of class \code{ExpressionSet}, longitudinal data object is created. Subsequantly the function estimates edges using partial correlation method with shrinkage approach applying \code{ggm.estimate.pcor} and \code{network.test.edges} functions. As a result an object or list of object of class \code{igraph} is created.
}

\usage{splineNetRecon(eSetObject, treatmentType, probesForNR="all", 
            cutoff.ggm=0.8, method="dynamic", saveEdges=FALSE)}

\arguments{
  \item{eSetObject}{\code{ExpressionSet} object of class \code{ExpressionSet} containing log-ratios or log-values of expression for a series of microarrays}
  \item{treatmentType}{a character string containing a type of \code{Treatment} defining samples considered for network reconstruction}
  \item{probesForNR}{a vector of character string containing names/IDs used for network reconstruction}
  \item{cutoff.ggm}{number or vector of numbers between 0 and 1 defining cutoff for significant posterior probability; default value is 0.8}
   \item{method}{method used to estimate the partial correlation matrix; available options are \code{"static"} and \code{"dynamic"} (default) - both are shrinkage methods}
   \item{saveEdges}{if \code{TRUE}, .Rdata file with all edges is created; default is \code{FALSE}}
}

\details{
The input \code{eSetObject} must be provided as an object of class \code{ExpressionSet} which contains \code{SampleName}, \code{Time}, \code{Treatment} and if applicable \code{Replicates} variables (columns) included in the phenotypic data of the \code{eSetObject} (\code{pData(eSetObject)}). Two types of \code{Treatment} defining two groups to compare have to be definied.

Gene association network reconstruction is conducted for a selected type of \code{Treatment}. This allows to find regulatory association between genes under a certain condition (treatment). First, a \code{longitudinal} data object of the gene expression data with possible repicates is created. This object is used to estimate partial correlation with selected shrinkage method (\code{"dynamic"} or \code{"static"}) with the \code{ggm.estimate.pcor} function (for details see \code{ggm.estimate.pcor} function help). Finally, the \code{network.test.edges} function estimates the probabilities for all possible edges and lists them in descending order (for details see \code{network.test.edges} help). 

\code{cutoff.ggm} can be a single number or a vector of numbers. If more than one value for \code{cutoff.ggm} is definied than function returns a \code{list} of objects of class \code{igraph} for each definied \code{cutoff.ggm} value. Otherwise a single object of class \code{igraph} with one selected probability is returned. 
}

\value{
An object or list of objects of class \code{igraph}.

If \code{saveEdges} is \code{TRUE}, .Rdata file with all possible edges is created.
}

\author{Agata Michna}
\seealso{
   \url{http://strimmerlab.org/software/genenet/}
   
   \url{http://strimmerlab.org/software/longitudinal/}
}

\examples{
## load "eSetObject" containing simulated time-course data
data(TCsimData)

## define function parameters 
treatmentType = "T2"
probesForNR = "all"
cutoff.ggm = 0.8
method = "dynamic"

## reconstruct gene association network from time-course data
igr <- splineNetRecon(eSetObject = TCsimData, treatmentType, probesForNR, cutoff.ggm, method)
plot(igr, vertex.label = NA, vertex.size = 3)
}

\keyword{differential expression}
\keyword{time-course data}
