#' bgamglm
#'
#' bgamglm code from: "Variational approximation for heteroscedastic models and matching
#' pursuit algorithms" by David Nott et al.
#'
#' @param w vector of responses
#' @param Z design matrix
#' @param pm prior mean
#' @param pv prior covariance matrix
#' @param p0 starting guess for mode
#' @param tol tolerance for stopping
#'
#' @author Minh-Ngoc Tran (ngoctm@nus.edu.sg) with a helper function from David Nott
#'
#' @export
#'
bgamglm <- function(w,Z,pm,pv,p0,tol) {
	# Helper function for the function vbhetlm1.
	# This function is to compute posterior mode for gamma generalized
	# linear model with log link, known shape parameters
	# and normal prior on coefficients
	#
	# w -
	# Z -
	# pm, pv -
	# p0 -
	# tolerance for stopping
	library(Matrix)
	xold = p0
	V = .5*w*exp(-Z%*%xold)
	fvalold = .5*sum(Z%*%xold)+sum(V)+0.5*sum((xold-pm)*solve(pv,(xold-pm)))
	Df = .5*colSums(Z)-t(Z)%*%V+solve(pv,xold-pm)
	W = Diagonal(x = as.numeric(V))
	D2f = t(Z)%*%W%*%Z+solve(pv)
	xnew = xold - solve(D2f,Df)
	V = .5*w*exp(-Z%*%xnew)
	fvalnew = .5*sum(Z%*%xnew)+sum(V)+0.5*sum((xnew-pm)*solve(pv,(xnew-pm)))
	fvalall = fvalold
	niter <- 1
	while((fvalold-fvalnew)>tol | niter < 20) {
		xold = xnew
		fvalold = fvalnew
		fvalall = c(fvalall,fvalnew)
		Df = .5*colSums(Z)-t(Z)%*%V+solve(pv,xold-pm)
		W = Diagonal(x = as.numeric(V))
		D2f = t(Z)%*%W%*%Z+solve(pv)
		xnew = xold - solve(D2f,Df)
		V = .5*w*exp(-Z%*%xnew)
		fvalnew = .5*sum(Z%*%xnew)+sum(V)+0.5*sum((xnew-pm)*solve(pv,(xnew-pm)))
		niter <- niter + 1
	}
	list(coef=xnew, variance=solve(D2f), fval=fvalall)
}
