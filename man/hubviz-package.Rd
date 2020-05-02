\name{hubviz-package}
\alias{hubviz-package}
\docType{package}

\title{
hubViz (A Novel Tool for Hub-centric Visualization)
}

\description{
This package provides functions for fitting hubViz,
 a novel tool for hub-centric visualization which is based on a latent space joint model (LSJM) for visualization. 
}

\details{
\tabular{ll}{
Package: \tab hubviz\cr
Type: \tab Package\cr
Version: \tab 0.1.0\cr
Date: \tab 2019-02-15\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}
This package contains a main class, \code{hubviz}, which represents hubViz model fit.

The main methods for the hubViz is \code{hubviz}, which fits a latent space joint model (LSJM).
}

\author{
Jin Hyun Nam, Ick Hoon Jin

Maintainer: Dongjun Chung <chungd@musc.edu>
}

\references{
Nam JH, Yun J, Jin IH and Chung D (2020) ''hubViz: A novel tool for hub-centric visualization''.
}

\keyword{ package }

\seealso{
\code{\link{hubviz}}
}


\examples{
library(hubviz)

# Data generation
n <- 100
p <- 10
bmat <- matrix( 0, n, p )
for ( j in 1:ncol(bmat) ) {
 bmat[ 1:(10*j), j ] <- 1
}
rownames( bmat ) <- paste( "Sample", 1:n, sep="" )
colnames( bmat ) <- paste( "V", 1:p, sep="" )

# Fit hubviz
nsamp <- nrow( bmat ); nitem <- ncol( bmat )
fit <- hubviz( bmat, nsample = nsamp, nitem = nitem)
fit

# Extraction of results
estimate( fit )

# Hub-centric plot
plot( fit )

# HPD plot
hpdplot( fit )
}