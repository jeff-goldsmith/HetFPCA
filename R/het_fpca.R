#' heteroskedastic FPCA
#'
#' @param X a design matrix for fixed effects in the mean model; it should have the same number of rows as Y. in an intercept-only
#' model this should be a column of ones, one for each curve
#' @param Y matrix of curves, one in each row
#' @param V design matrix for random effects; it should have the same number of rows as Y. a 1 row i, column j indicates that curve
#' i is made by subject j. use V=NULL for no random effects in the mean model
#' @param het.model.matrix design matrix for the score variance model; it should have the same number of rows as Y. each column
#' corresponds to a score variance covariate. use het.model.matrix=NULL for no heteroskedasticity
#' @param het.model.prior.component with length ncol(het.model.matrix). each element indicates to which variance component the
#' corresponding covariate belongs.
#' @param fix.sigma used to indicate that a covariate in the score variance model should be a fixed effect, not a random effect. it
#' should be a data frame with two columns, comp and val: comp corresponds to a component in het.model.prior.component; val is the
#' value of the normal prior, e.g., 100, to assign to that component
#' @param Kt number of spline basis functions
#' @param Kp number of FPCs to estimate
#' @param ig.alpha parameter for the inverse gamma priors used to control smoothness of functions and for the variance of the
#' random effects in the score variance model
#' @param ig.beta parameter for the inverse gamma priors used to control smoothness of functions and for the variance of the
#' random effects in the score variance model
#' @param max.iterations max number of iterations
#' @param tol tolerance
#'
#' @return
#'
#' @export
#'
het_fpca <- function (X, Y=NULL, V=NULL,
          het.model.matrix=NULL, het.model.prior.component=NULL,
          fix.sigma=NULL,
          Kt = 5, Kp = 2,
          ig.alpha=1, ig.beta=1,
          max.iterations=10, tol=5e-6)
{
  I = dim(X)[1]         # number of curves
  D = dim(Y)[2]         # length of curves
  p = dim(X)[2]         # number of predictors
  if (is.null(V)){      # if there are no random effects make dummy matrix
    ranefs <- F
    V <- matrix(0, nrow=I, ncol=1)
  } else {
    ranefs <- T
  }
  heteroskedastic <- !is.null(het.model.matrix)
  library(splines)
  Theta = bs(1:D, df = Kt, intercept = TRUE, degree = 3)
  diff2 = matrix(rep(c(1, -2, 1, rep(0, D - 2)),
                     D - 2)[1:((D - 2) * D)],
                 D - 2, D, byrow = TRUE)
  P.mat = t(Theta) %*% t(diff2) %*% diff2 %*% Theta
  Y.vec = as.vector(t(Y))
  obspts.vec = !is.na(Y.vec)
  Y.vec = Y.vec[obspts.vec]
  J = sum(obspts.vec)
  t.designmat.X = t(kronecker(X, Theta)[obspts.vec, ])
  XtX = matrix(0, Kt * p, Kt * p)
  sumXtX = matrix(0, Kt * p, Kt * p)
  for (i in 1:I) {
    obs.points = which(!is.na(Y[i, ]))
    X.cur = kronecker(matrix(X[i, ], nrow = 1, ncol = p),
                      Theta)[obs.points, ]
    XtX = XtX + crossprod(X.cur)
    sumXtX = sumXtX + t(X.cur) %*% X.cur
  }
  # initialization
  mu.q.BW <- matrix(0, nrow=Kt, ncol=p)
  b.q.lambda.BW = rep(ig.beta, p)
  mu.q.BZ <- matrix(0, nrow=Kt, ncol=ncol(V))
  b.q.lambda.BZ <- ig.beta
  sigma.q.Bpsi = vector("list", Kp)
  for (k in 1:Kp) {
    sigma.q.Bpsi[[k]] = diag(1, Kt)
  }
  mu.q.Bpsi = matrix(0, nrow = Kt, ncol = Kp)
  b.q.lambda.Bpsi = rep(ig.beta, Kp)
  sigma.q.C = vector("list", I)
  for (k in 1:I) {
    sigma.q.C[[k]] = diag(1, Kp)
  }
  set.seed(2)    # for reproducibility of results, since algorithm
                 # initialized with random scores
  mu.q.C = matrix(rnorm(I * Kp, 0, 0.01), I, Kp)
  b.q.sigma.me = ig.beta
  pcaef.cur = matrix(0, I, D)
  j <- 1
  inv.sigma.me <- (ig.alpha + J /2)/(b.q.sigma.me)
  while (j < max.iterations) {
    # estimate mean
    mean.cur = as.vector(t(pcaef.cur) + t(V %*% t(mu.q.BZ) %*% t(Theta)))[obspts.vec]
    sigma.q.beta = solve(inv.sigma.me *
                           XtX +
                           kronecker(diag((ig.alpha + Kt/2)/b.q.lambda.BW,
                                          length(b.q.lambda.BW)), P.mat))
    mu.q.beta = matrix(sigma.q.beta %*%
                         (inv.sigma.me * t.designmat.X %*% (Y.vec - mean.cur)),
                       nrow = Kt, ncol = p)
    fixef.cur = as.matrix(X %*% t(mu.q.beta) %*% t(Theta))
    # estimate random effects in mean model
    if (ranefs){
      exp.bT.Q.b <- 0
      for (i in 1:ncol(V)){
        which.curves <- which(V[, i]==1); num.curves <- length(which.curves)
        theta.mats <- vector('list', num.curves)
        for (curve in 1:num.curves){
          theta.mats[[curve]] <- Theta[!is.na(Y[which.curves[curve], ]), ]
        }
        design.mat <- do.call('rbind', theta.mats)
        sigma.q.BZ <- solve(inv.sigma.me *
                              t(design.mat) %*% design.mat +
                              kronecker(diag((ig.alpha + ncol(V) * Kt/2)/b.q.lambda.BZ,
                                             length(b.q.lambda.BZ)), P.mat))
        Y.subj <- Y[which.curves, ]
        obspts.vec <- !is.na(as.vector(t(Y.subj)))
        mean.cur = as.vector(t(fixef.cur[which.curves, ] +
                               pcaef.cur[which.curves, ]))[obspts.vec]
        Y.subj.vec <- as.vector(t(Y.subj))[obspts.vec]
        mu.q.BZ[, i] = sigma.q.BZ %*% (inv.sigma.me * t(design.mat)) %*%
                                            (Y.subj.vec-mean.cur)
        for (term in 1){
          indices <- (nrow(P.mat) * (term - 1) + 1):(nrow(P.mat) * term)
          exp.bT.Q.b[term] <- exp.bT.Q.b[term] + t(mu.q.BZ[indices, i]) %*% P.mat %*% mu.q.BZ[indices, i] +
            sum(diag(P.mat %*% sigma.q.BZ[indices, indices]))
        }
      }
    }
    ranef.cur = V %*% t(mu.q.BZ) %*% t(Theta)
    designmat = kronecker(mu.q.C, Theta)[obspts.vec, ]
    # estimate FPCs
    sigma.q.Bpsi = solve(kronecker(diag((ig.alpha + Kt/2)/b.q.lambda.Bpsi),
                                   P.mat) + inv.sigma.me *
                           f_sum(mu.q.c = mu.q.C, sig.q.c = sigma.q.C, theta = t(Theta),
                                 obspts.mat = !is.na(Y)))
    mu.q.Bpsi = matrix(inv.sigma.me * sigma.q.Bpsi %*%
                         f_sum2(y = Y, fixef = fixef.cur + ranef.cur, mu.q.c = mu.q.C,
                                kt = Kt, theta = t(Theta)), nrow = Kt, ncol = Kp)
    # orthogonalize FPCs
    psi.cur <- svd(mu.q.C %*% t(mu.q.Bpsi) %*% t(Theta))$v[, 1:Kp]
    psi.cur <- apply(psi.cur, 2, function(x){x/sign(x[1])})
    psi.cur <- t(psi.cur)
    # estimate scores
    if (j==1 & heteroskedastic){
      inverse.score.variances <- matrix(1, I, Kp)
    }
    for (subj in 1:I) {
      obs.points = which(!is.na(Y[subj, ]))
      Theta_i = t(Theta)[, obs.points]
      if (is.null(het.model.matrix)){
        m <- diag(1, Kp, Kp)
      } else {
        m <- diag(inverse.score.variances[subj, ], Kp, Kp)
      }
      sigma.q.C[[subj]] = solve(m + inv.sigma.me * (f_trace(Theta_i = Theta_i,
                      Sig_q_Bpsi = sigma.q.Bpsi, Kp = Kp, Kt = Kt) +
                        psi.cur[, obs.points] %*% t(psi.cur[, obs.points])))
      mu.q.C[subj, ] = inv.sigma.me *
        sigma.q.C[[subj]] %*% as.matrix(psi.cur[, obs.points]) %*%
        (Y[subj, obs.points] - fixef.cur[subj, obs.points] - ranef.cur[subj, obs.points])
    }
    # initialize score variance model parameters
    if (!heteroskedastic){
      g_mean_l <- g_cov_l <- NULL
    } else {
      var.p <- ncol(het.model.matrix)
      if (j==1) {
        reg.inverse.score.variances <- rep(1, Kp)
        prior_cov_l <- g_mean_l <- g_cov_l <- p0 <- vector('list', Kp)
        for (k in 1:Kp){
          g_mean_l[[k]] <- rep(1, var.p)
          g_cov_l[[k]] <- diag(1, var.p)
          prior_cov_l[[k]] <- diag(1, var.p)
          p0[[k]] <- rep(0, var.p)
        }
      } else {
        p0 <- g_mean_l
      }
      # estimate parameters in score variance model
      for (k in 1:Kp){
        first.order.term <- mu.q.C[, k]^2
        second.order.term <- unlist(lapply(sigma.q.C, function(x){x[k,k]}))
        mu.C_squared <- first.order.term + second.order.term
        gr <- GammaRegression(mu.C_squared=mu.C_squared,
                              prior.cov=prior_cov_l[[k]],
                              p0=p0[[k]],
                              het.model.matrix=het.model.matrix,
                              het.model.prior.component=het.model.prior.component,
                              fix.sigma=fix.sigma,
                              ig.alpha=1, ig.beta=1)
        g_mean_l[[k]] <- gr$g_mean
        g_cov_l[[k]] <- gr$g_cov
        prior_cov_l[[k]] <- gr$het.model.prior.var
        inverse.score.variances[, k] <- gr$recip.exp.g.x
      }
    }
    pcaef.cur = as.matrix(mu.q.C %*% psi.cur)
    resid = as.vector(Y - fixef.cur - pcaef.cur - ranef.cur)
    # estimate variances
    b.q.sigma.me = as.numeric(ig.beta +
                                0.5 * (crossprod(resid[!is.na(resid)]) +
                                         sum(diag(sumXtX %*% sigma.q.beta)) +
                                         f_sum4(mu.q.c = mu.q.C,
                                                sig.q.c = sigma.q.C,
                                                mu.q.bpsi = mu.q.Bpsi,
                                                sig.q.bpsi = sigma.q.Bpsi,
                                                theta = Theta,
                                                obspts.mat = !is.na(Y),
                                                psi.cur=psi.cur)))
    inv.sigma.me <- (ig.alpha + J /2)/(b.q.sigma.me)
    for (term in 1:length(b.q.lambda.BW)) {
      mu.indices <- 1:(nrow(P.mat))
      s.indices <- (nrow(P.mat) * (term - 1) + 1):(nrow(P.mat) * term)
      col.num <- term
      b.q.lambda.BW[term] = ig.beta + 0.5 *
        (t(mu.q.BW[mu.indices, col.num]) %*% P.mat %*%
          mu.q.BW[mu.indices, col.num] +
           sum(diag(P.mat %*% sigma.q.beta[s.indices, s.indices])))
    }
    for (term in 1:length(b.q.lambda.Bpsi)) {
      mu.indices <- 1:nrow(P.mat)
      s.indices <- (nrow(P.mat) * (term - 1) + 1):(nrow(P.mat) * term)
      col.num <- term
      b.q.lambda.Bpsi[term] = ig.beta +
        0.5 * (t(mu.q.Bpsi[mu.indices, col.num]) %*% P.mat %*% mu.q.Bpsi[mu.indices, col.num] +
        sum(diag(P.mat %*% sigma.q.Bpsi[s.indices, s.indices])))
    }
    if (ranefs){
      for (term in 1){
        b.q.lambda.BZ[term] = as.numeric(ig.beta + 0.5 * exp.bT.Q.b[term])
      }
    }
    # monitor convergence
    current.objective.mean <- c(as.numeric(mu.q.beta, mu.q.Bpsi, mu.q.C))
    current.objective.var <- unlist(g_mean_l)
    if (j >= 2){
      diff.mean <- sum((current.objective.mean-old.objective.mean)^2)
      tol.mean <- diff.mean < tol
      if (!is.null(het.model.matrix)){
        diff.var <- sum((current.objective.var-old.objective.var)^2)
        tol.var <- diff.var < tol
      }
      if (is.null(het.model.matrix)){
        if (tol.mean) break
        cat("Diff.mean", diff.mean, "\n")
      } else {
        if (tol.mean & tol.var & j > 5) break
        cat("Diff.mean", diff.mean, "Diff.var", diff.var, "\n")
      }
    }
    old.objective.mean <- current.objective.mean
    old.objective.var <- current.objective.var
    j <- j + 1
  }
  if (!is.null(het.model.matrix)){
    g_mean <- do.call('cbind', g_mean_l)
  } else {
  	g_mean = NULL
  }
  ret = list(fixef = fixef.cur,
  					 ranef = ranef.cur,
  					 pcaef = pcaef.cur,
  					 psi.cur=psi.cur,
             g_mean=g_mean, g_cov=g_cov_l,
             mu.q.beta=mu.q.beta, mu.q.Bpsi=mu.q.Bpsi, mu.q.BZ=mu.q.BZ,
             sigma.q.beta=sigma.q.beta, sigma.q.Bpsi=sigma.q.Bpsi,
             mu.q.C=mu.q.C, Theta=Theta, het.model.matrix=het.model.matrix,
             het.model.prior.component=het.model.prior.component, Kp=Kp, Kt=Kt,
             heteroskedastic=heteroskedastic)
  ret
}

f_sum = function(mu.q.c, sig.q.c, theta, obspts.mat){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  kt = dim(theta)[1]
  ret.sum = matrix(0, kp*kt, kp*kt)
  for(i in 1:I){
    mu.mat = matrix(mu.q.c[i,], nrow = 1, ncol = kp)
    ret.sum = ret.sum + kronecker(t(mu.mat) %*% mu.mat + sig.q.c[[i]], (theta[,obspts.mat[i,]])%*%t(theta[,obspts.mat[i,]]))
  }
  return(ret.sum)
}

f_sum2 = function(y, fixef, mu.q.c, kt, theta){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  ret.sum = matrix(0, nrow = kp*kt, ncol = 1)
  for(i in 1:I){
    obs.pts = !is.na(y[i,])
    ret.sum = ret.sum + kronecker((matrix(mu.q.c[i,])), theta[,obs.pts]) %*% matrix(y[i, obs.pts] - fixef[i,obs.pts])
  }
  return(ret.sum)
}

f_sum4 = function(mu.q.c, sig.q.c, mu.q.bpsi, sig.q.bpsi, theta, obspts.mat, psi.cur){
  I = dim(mu.q.c)[1]
  kp = dim(mu.q.c)[2]
  kt = dim(theta)[2]
  ret.sum = matrix(0, 1, 1)
  for(i in 1:I){
    theta_i = t(theta)[,obspts.mat[i,]]
    t.mu.q.bpsi_theta_i <- psi.cur[, obspts.mat[i, ]]
    temp =
      f_trace(Theta_i = theta_i, Sig_q_Bpsi = sig.q.bpsi, Kp = kp, Kt = kt) %*% matrix(mu.q.c[i,], kp, 1) %*% matrix(mu.q.c[i,], 1, kp) +
      f_trace(Theta_i = theta_i, Sig_q_Bpsi = sig.q.bpsi, Kp = kp, Kt = kt) %*% sig.q.c[[i]] +
      t.mu.q.bpsi_theta_i %*% t(t.mu.q.bpsi_theta_i) %*% sig.q.c[[i]]
    ret.sum = ret.sum + sum(diag(temp))
  }
  return(ret.sum)
}

f_trace = function(Theta_i, Sig_q_Bpsi, Kp, Kt){
  ret.mat = matrix(NA, nrow = Kp, ncol = Kp)
  A = Theta_i %*% t(Theta_i)
  for(i in 1:Kp){
    for(j in 1:Kp){
      ret.mat[i,j] = sum(diag(A %*% Sig_q_Bpsi[((-1 + i)*Kt + 1):(i*Kt), ((-1 + j)*Kt + 1):(j*Kt)]))
    }
  }
  return(ret.mat)
}

GammaRegression <- function(mu.C_squared, prior.cov=100, p0,
                            het.model.matrix,
                            het.model.prior.component=NULL,
                            fix.sigma=NULL,
                            ig.alpha=1, ig.beta=1){
  mu.C_squared <- matrix(mu.C_squared)
  prior.mean <- matrix(rep(0, ncol(het.model.matrix)))
  fit <- bgamglm(mu.C_squared, het.model.matrix,
                      prior.mean, prior.cov,
                      p0=p0, tol=1e-6)
  fit.recip.exp.g.x <- exp(
    -het.model.matrix %*% fit$coef
    -1/2 * diag(het.model.matrix %*%
                  fit$variance %*%
                  t(het.model.matrix) ))
  g_mean <- as.numeric(fit$coef)
  recip.exp.g.x <- as.numeric(fit.recip.exp.g.x)
  g_cov <- fit$variance
  if (!is.null(het.model.prior.component)){
    het.model.prior.var <-
      matrix(0, nrow=length(het.model.prior.component),
             ncol=length(het.model.prior.component))
    num.fixed <- nrow(fix.sigma)
    for (l in 1:num.fixed){
      which.indices <- which(het.model.prior.component==
                               fix.sigma[l, "comp"])
      if (length(which.indices)>1){
        diag(het.model.prior.var[
          which.indices, which.indices]) <-
          fix.sigma[l, "val"]
      } else {
        het.model.prior.var[which.indices, which.indices] <-
          fix.sigma[l, "val"]
      }
    }
    random.components <- setdiff(unique(het.model.prior.component),
                                 fix.sigma$comp)
    for (k in 1:length(random.components)){
      w1 <- which(het.model.prior.component==random.components[k])
      s <- UpdateGSigma(
        indices=w1,
        g_mean=g_mean, g_cov=g_cov,
        ig.alpha=ig.alpha, ig.beta=ig.beta)
      het.model.prior.var[cbind(w1,w1)] <- s
    }
  } else {
    het.model.prior.var=NULL
  }
  return(list(g_mean=g_mean, g_cov=g_cov, recip.exp.g.x=recip.exp.g.x,
              het.model.prior.var=het.model.prior.var))
}
