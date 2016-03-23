\name{splineDiffExprs}
\alias{splineDiffExprs}
\title{Differential expression analysis based on natural cubic spline regression models for time-course data}

\description{
The function compares time dependent behaviour of genes in two different groups. Applying empirical Bayes moderate F-statistic on differences in coefficients of fitted natural cubic spline regression models, differentially expressed in time genes are determined. The function is a wrapper of other R-functions to simplify differential expression analysis of time-course data.
}

\usage{splineDiffExprs(eSetObject, df, cutoff.adj.pVal=1, reference, intercept=TRUE)}

\arguments{
  \item{eSetObject}{\code{ExpressionSet} object of class \code{ExpressionSet} containing log-ratios or log-values of expression for a series of microarrays}
  \item{df}{number of degrees of freedom}
  \item{cutoff.adj.pVal}{Benjamini-Hochberg adjusted p-value cut-off}
  \item{reference}{character defining which treatment group should be considered as reference}
  \item{intercept}{if \code{TRUE}, F-test includes all parameters; if \code{FALSE}, F-test includes shape parameters only; default is \code{TRUE}}
}

\details{
The function fits a temporal trend using a natural cubic spline regression to simulate nonlinear behaviour of genes over time. 

The input \code{eSetObject} must be provided as an object of class \code{ExpressionSet} which contains \code{SampleName}, \code{Time}, \code{Treatment} and if applicable \code{Replicates} variables (columns) included in the phenotypic data of the \code{eSetObject} (\code{pData(eSetObject)}). Two types of \code{Treatment} defining two groups to compare have to be definied. 

Replicates are not required. The time points for compared treatment groups should be identical.

User has to define number of degrees of freedom (\code{df}) for the spline regression model. Choosing effective degrees of freedom in range 3-5 is reasonable. 

Time dependent differential expression of a gene is determined by the application of empirical Bayes moderate F-statistics on the differences of coefficient values of the fitted natural cubic spline regression models for the same gene in the two compared treatment groups. In other words, comparing the coefficient values of the fitted splines in both groups allows the detection of differences in the shape of the curves, which represent the gene expressions changes over time. Ouptut table containing Benjamini-Hochberg adjusted p-value (\code{adj.P.Value}) is used to define differentially expressed genes. The default value for \code{cutoff.adj.pVal} is set to \code{1}, which means that all genes are included in output table.
}
\value{
A data.frame with rows defining names/IDs of differentially expressed genes and additional columns described below.

The first columns contain all feature data of the \code{eSetObject} (\code{fData(eSetObject)}), if any feature data were defined. Otherwise, only one column \code{row_IDs}, containing the row names is created. The \code{b_0}, \code{b_1},..., \code{b_m} coefficients correspond to the reference model parameters. The \code{d_0}, \code{d_1},..., \code{d_m} coefficients represent the differences between the reference model parameters and the model parameters in the compared group. \code{AveExprs} refers to the average log2-expression for a probe (representing a gene) over all arrays. The \code{F} column contains moderate F-statistics, \code{P.Value} raw p-value and \code{adj.P.Value} Benjamini-Hochberg adjusted p-value.
}

\author{Agata Michna}
\seealso{\code{\link{limma}}}

\examples{
## load "eSetObject" containing simulated time-course data
data(TCsimData)
pData(TCsimData)

## define function parameters 
df <- 3
cutoff.adj.pVal <- 0.01
reference <- "T1"
intercept <- TRUE

diffExprs <- splineDiffExprs(eSetObject = TCsimData, df, cutoff.adj.pVal, reference, intercept)
head(diffExprs,3)
}

\keyword{differential expression}
\keyword{time-course data}
\keyword{spline}

