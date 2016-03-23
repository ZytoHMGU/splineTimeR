\name{splinePlot}
\alias{splinePlot}
\title{Plot spline regression curves of time-course data}

\description{
Function visualises time dependent behaviour of genes in two compared groups. The natural cubic spline regression curves fitted to discrete, time dependent expression data are plotted. One plot shows two curves - representing the reference group and the compared group, respectively. See also \code{splineDiffExprs} function.
}

\usage{splinePlot(eSetObject, df, reference, toPlot="all")}

\arguments{
  \item{eSetObject}{\code{ExpressionSet} object of class \code{ExpressionSet} containing log-ratios or log-values of expression for a series of microarrays}
  \item{df}{number of degrees of freedom}
  \item{reference}{character defining which treatment group should be considered as reference}
  \item{toPlot}{vector of genes to plot; defalut is \code{toPlot} = \code{"all"}}
}

\details{
The input \code{eSetObject} must be provided as an object of class \code{ExpressionSet} which contains \code{SampleName}, \code{Time}, \code{Treatment} and if applicable \code{Replicates} variables (columns) included in the phenotypic data of the \code{eSetObject} (\code{pData(eSetObject)}). Two types of \code{Treatment} defining two groups to compare have to be definied.

Replicates are not required. The time points for compared treatment groups should be identical.

User has to define number of degrees of freedom (\code{df}) for the spline regression model. Choosing effective degrees of freedom in range 3-5 is reasonable. 

Genes to plot, given as a vector of characters, can be selected by the user. Provided names have to be a part of a row name vector of \code{eSetObject} (\code{rownames(exprs(eSetObject))}). If genes to plot are not definied, all genes are plotted. 
}

\value{A .pdf file containing plots for chosen genes.}
\author{Agata Michna}
\seealso{\code{\link{limma}}}

\examples{
## load "eSetObject" object containing simulated time-course data
data(TCsimData)
pData(TCsimData)

## define function parameters 
df <- 3
reference <- "T1"
toPlot <- rownames(TCsimData)[1:10]

splinePlot(eSetObject = TCsimData, df, reference, toPlot)
}

\keyword{spline}
\keyword{time-course data}
