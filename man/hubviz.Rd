
\name{hubviz}
\alias{hubviz}

\title{hubViz: A Novel Tool for Hub-centric Visualization}

\description{
Fit hubViz model
}

\usage{
hubviz(dataset, nsample = nrow(dataset), nitem = ncol(dataset),
  ndim = 2, niter = 30000, nburn = 5000, nthin = 5, nprint = 100,
  jump_theta = 1, jump_w = 0.2, pr_mean_theta = 0,
  pr_sd_theta = 10, pr_mean_w = 0, prior_a = 0.001,
  prior_b = 0.001, w.ci=0.5, option = FALSE, verbose = FALSE, cores = 1)
}

\arguments{
\item{dataset}{Matrix of binary values with n samples and p variables}

\item{nsample}{Number of samples to analyze (default: use all samples)}

\item{nitem}{Number of variables to analyze (default: use all variables)}

\item{ndim}{Dimension of latent position (default: 2)}

\item{niter}{Total number of iterations in MCMC (default: 30000)}

\item{nburn}{Number of burn-in iterations in MCMC (default: 5000)}

\item{nthin}{Thinning number of MCMC (default: 5)}

\item{nprint}{Trace size of MCMC (default: 100)}

\item{jump_theta}{Jump size for the intercept parameter (default: 1.0)}

\item{jump_w}{Jump size for latent space (default: 0.2)}

\item{pr_mean_theta}{Prior mean for the intercept parameter (default: 0.0)}

\item{pr_sd_theta}{Prior Standard deviation for the intercept parameter (default: 10.0)}

\item{pr_mean_w}{Prior mean for latent space (default: 0.0)}

\item{prior_a}{Prior shape parameter for variance of latent space (default: 0.001)}

\item{prior_b}{Prior scale parameter for variance of latent space (default: 0.001)}

\item{w.ci}{Create highest posterior density (HPD) intervals for latent positions in an MCMC (default: 0.5)}

\item{option}{Logical. If TRUE then print MCMC iteration (default: FALSE)}

\item{verbose}{Logical. If TRUE then return MCMC results (default: FALSE)}

\item{cores}{Number of CPUs to be used for computing (default: 1)}
}

\details{
\code{hubviz} fits the hubViz model and requires binary matrix as input data.
}

\value{
Construct \code{hubviz} class object.
}

\references{
Nam JH, Yun J, Jin IH and Chung D (2020) ''hubViz: A novel tool for hub-centric visualization''.
}

\author{ Jin Hyun Nam, Ick Hoon Jin }


\examples{
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

# Extract latent positions
estimate( fit )

# Hub-centric plot
plot( fit )

# HPD plot
hpdplot( fit )
}
