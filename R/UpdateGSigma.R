#' UpdateGSigma
#'
#' @param indices
#' @param g_mean
#' @param g_cov
#' @param ig.alpha
#' @param ig.beta
#'
#' @export
#'
UpdateGSigma <- function(indices, g_mean, g_cov, ig.alpha, ig.beta){
	alpha <- length(indices)/2 + ig.alpha
	coef <- g_mean[indices]
	trace.cov <- sum(diag(g_cov[indices, indices, drop = FALSE]))
	beta <- 1/2 * (sum(coef * coef) + trace.cov) + ig.beta
	beta / alpha
}
