\name{networkProperties}
\alias{networkProperties}
\title{Scale-free properties of a network}

\description{
Function plots degree distribution of nodes in a given network.
}
\usage{networkProperties(igr)}
\arguments{
  \item{igr}{an object or list of objects of class \code{igraph}}
}

\details{
Biological networks are thought to be scale-free. Scale-free networks follow a power-law distribution of the degrees of nodes in the network. This distribution is characterised by the degree exponent gamma, which for biological networks ranges between 2 and 3.

The function calculates the degree exponent(s) of given network(s) in comparison with degree exponents of biological networks derived from Reactome and BioGRID repositories. Both of these networks are built using functional interaction pairs extracted from mentioned repositories and provided in \code{FIs} data package (see example).

For each given \code{igraph} three types of plots are created: empirical cumulative distribution, degree distribution and power-law degree distribution on log-log scale with fitted trend line. 
}

\value{
A summary table containing number of nodes, number of edges and degree exponents for each given network.

Additionally two .pdf files are created. One containing empirical cumulative distribution frequency plots together with degree distributions and second with plots of power-law degree distribution on log-log scale.
}

\author{Agata Michna, Martin Selmansberger}
\seealso{\code{\link{FIs}}}

\references{
Barabasi, A. L. and Albert, R. (1999). Emergence of scaling in random networks.
\emph{Science} 286, 509-512.

Barabasi, A. L. and Oltvai, Z. N. (2004). Network biology: understanding the cell's functional organization.
\emph{Nat. Rev. Genet.} 5, 101-113.
}

\examples{
## load "eSetObject" containing simulated time-course data
data(TCsimData)

## reconstruct gene association networks from time-course data
igr <- splineNetRecon(eSet = TCsimData, treatmentType = "T2", cutoff.ggm = c(0.8,0.9))

## check for scale-free properties of reconstructed networks (igraphs)
scaleFreeProp <- networkProperties(igr)
head(scaleFreeProp)

## the functional interaction pairs provided in FIs data package
library(FIs)
data(FIs)
names(FIs)
head(FIs$FIs_Reactome)
head(FIs$FIs_BioGRID)
}

\keyword{scale-free properties}
\keyword{degree distribution}
