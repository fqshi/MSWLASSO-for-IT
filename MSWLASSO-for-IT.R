#' Solve cardinality contrained index tracking problem
#'
#' @param y The index return.
#' @param R The matrix of stock returns.
#' @param k The number of selected stocks.
#' @param l The lower bound of weight.
#' @param u The upper bound of weight.
#' @param TimeLimit Time limit.
#' @param LogFile
#' @param OutputFlag
#' @return The optimal portfolio weights.
card_gurobi <- function(y, R, k, l = 0, u = 1, TimeLimit = 1200) {
  p <- ncol(R)
  if (length(l) == 1) {
    l <- rep(l, p)
  }
  if (length(u) == 1) {
    u <- rep(u, p)
  }
  Qmat <- rbind(cbind(t(R) %*% R, matrix(0, ncol = p, nrow = p)),
    matrix(0, nrow = p, ncol = 2 * p))
  cvec <- c(-2 * as.vector(y %*% R), rep(0, p))
  Amat_up <- cbind(diag(rep(1, p)), diag(-u))
  Amat_mid <- cbind(diag(rep(-1, p)), diag(l))
  Amat_down <- rbind(cbind(matrix(1, nrow = 1, ncol = p),
                           matrix(0, nrow = 1, ncol = p)),
    cbind(matrix(0, nrow = 1, ncol = p), matrix(1, nrow = 1, ncol = p)))
  Amat <- rbind(Amat_up, Amat_mid, Amat_down)
  bvec <- c(rep(0, 2 * p), 1, k)
  model <- list()
  model$Q <- Qmat
  model$obj <- cvec
  model$A <- Amat
  model$rhs <- bvec
  model$sense <- c(rep("<", 2 * p), "=", "=")
  model$vtypes <- c(rep("C", p), rep("B", p))
  model$start <- c(rep(1 / p, p), rep(NA, p))

  params <- list()
  params$TimeLimit <- TimeLimit
  params$LogToConsole = 0
  params$OutputFlag = 0
  res <- gurobi(model, params)
  w <- res$x
  return(w[1: (length(w) / 2)])
}

#' One-step LLA
#'
#' @param y Index returns.
#' @param x Stock returns.
#' @param lambda Tunning parameter.
#' @param aw_fun Weighting function.
#' @param a Parameter used in aw_fun.
#' @param ub Upper bound.
#' @param tolerance
#' @return
alit_one <- function(y, x, lambda, aw_fun, a, ub = 1, tol = 1e-6) {
  p <- ncol(x)
  w0 <- qpit(y, x, aw = rep(1, p), ub = ub, tol = tol)
  aw_fun <- match.fun(aw_fun)
  aw <- aw_fun(w0, lambda, a)
  index_finite <- which(is.finite(aw))
  x_subset <- x[, index_finite, drop = FALSE]
  aw_subset <- aw[index_finite]
  w_subset <- qpit(y, x_subset, aw = aw_subset, ub = ub, tol = tol)
  w_whole <- rep(0, p)
  w_whole[index_finite] <- w_subset
  return(w_whole)
}


#' Multi-step LLA
alit_multi <- function(y, x, lambda, aw_fun, a, ub = 1, tol = 1e-6) {
  p <- ncol(x)
  w_iter <- qpit(y, x, aw = rep(1, p), ub = ub, tol = tol)
  aw_fun <- match.fun(aw_fun)
  niter <- 0
  while (TRUE) {
    if (niter > 500) {
      print("alit_multi: iteration times > 500")
      return(w_iter)
    } else {
      niter <- niter + 1
    }
    aw <- aw_fun(w_iter, lambda, a)
    index_finite <- which(is.finite(aw))
    x_subset <- x[, index_finite, drop = FALSE]
    aw_subset <- aw[index_finite]
    w_subset <- qpit(y, x_subset, aw = aw_subset, ub = ub, tol = tol)
    w_whole <- rep(0, p)
    w_whole[index_finite] <- w_subset
    if (all(abs(w_whole - w_iter) < tol)) {
      w <- w_whole
      break
    } else {
      w_iter <- w_whole
    }
  }
  return(w)
}


#' Quadratic programming for index tracking with adptive lasso
#'
#' @param y A vector of index return.
#' @param x A n \times p matrix of asset returns.
#' @param aw A vector of weights.
#' @param ub A positive scalar upper bound. Default 1.
#' @param tol Default 1e-6.
#'
#' @return A vector of the optimal weights.
#'
#' @export
qpit <- function(y, x, aw, ub = 1, tol = 1e-6) {
  p <- ncol(x)
  n <- nrow(x)
  Vmat <- t(x) %*% x * 2
  dvec <- aw - as.vector(t(y) %*% x) * 2
  Amat <- matrix(1, nrow = 1, ncol = p)
  bvec <- 1
  uvec <- rep(ub, p)
  method <- "CHOL"
  f <- file()
  sink(file = f)
  res <- LowRankQP::LowRankQP(Vmat = Vmat, dvec = dvec, Amat = Amat,
                              bvec = bvec, uvec = uvec, method = method)
  # for silencing message
  sink()
  close(f)
  w <- as.vector(res$alpha)
  w[w < tol] <- 0
  return(w)
}

aw_scad <- function(u, lambda, a) {
  if (a <= 2) {
    stop("a should be larger than 2.")
  }
  u_abs <- abs(u)
  pen_prime_abs <- (lambda * (u_abs <= lambda)
    + pmax(a * lambda - u_abs, 0) / (a - 1) * (u_abs > lambda))
  return(pen_prime_abs)
}

aw_mcp <- function(u, lambda, a) {
  if (a <= 1) {
    stop("a should be larger than 1.")
  }
  pen_prime_abs <- pmax((a * lambda - abs(u)) / a, 0)
  return(pen_prime_abs)
}

aw_lq <- function(u, lambda, a) {
  if (a >= 1 || a <= 0) {
    stop("a should be in (0, 1)")
  }
  pen_prime_abs <- ((abs(u))^(a - 1)) * a * lambda
  return(pen_prime_abs)
}

aw_log_m <- function(u, lambda, a) {
  if (a <= 0) {
    stop("s should be larger than 0.")
  }
  pen_prime_abs <- lambda / log(1 + 1 / a) / (a + abs(u))
  return(pen_prime_abs)
}
