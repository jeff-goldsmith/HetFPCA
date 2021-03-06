% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/het_fpca.R
\name{het_fpca}
\alias{het_fpca}
\title{heteroskedastic FPCA}
\usage{
het_fpca(X, Y = NULL, V = NULL, het.model.matrix = NULL,
  het.model.prior.component = NULL, fix.sigma = NULL, Kt = 5, Kp = 2,
  ig.alpha = 1, ig.beta = 1, max.iterations = 10, tol = 5e-06)
}
\arguments{
\item{X}{a design matrix for fixed effects in the mean model; it should have the same number of rows as Y. in an intercept-only
model this should be a column of ones, one for each curve}

\item{Y}{matrix of curves, one in each row}

\item{V}{design matrix for random effects; it should have the same number of rows as Y. a 1 row i, column j indicates that curve
i is made by subject j. use V=NULL for no random effects in the mean model}

\item{het.model.matrix}{design matrix for the score variance model; it should have the same number of rows as Y. each column
corresponds to a score variance covariate. use het.model.matrix=NULL for no heteroskedasticity}

\item{het.model.prior.component}{with length ncol(het.model.matrix). each element indicates to which variance component the
corresponding covariate belongs.}

\item{fix.sigma}{used to indicate that a covariate in the score variance model should be a fixed effect, not a random effect. it
should be a data frame with two columns, comp and val: comp corresponds to a component in het.model.prior.component; val is the
value of the normal prior, e.g., 100, to assign to that component}

\item{Kt}{number of spline basis functions}

\item{Kp}{number of FPCs to estimate}

\item{ig.alpha}{parameter for the inverse gamma priors used to control smoothness of functions and for the variance of the
random effects in the score variance model}

\item{ig.beta}{parameter for the inverse gamma priors used to control smoothness of functions and for the variance of the
random effects in the score variance model}

\item{max.iterations}{max number of iterations}

\item{tol}{tolerance}
}
\value{

}
\description{
heteroskedastic FPCA
}
