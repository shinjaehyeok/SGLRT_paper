#' Upper bound on the boundary crossing probability for confidence sequences.
#'
#' \code{cross_prop_cs} is used to compute the constant threshold for the GLR-like confidence sequences with finite target time interval.
#'
#' @param g Constant threshold (positive numeric).
#' @param nmax Upper bound of the target time interval
#' @param nmin Lower bound of the target time interval (default = 1L)
#' @param eta Peeling size (default = 1.1).
#'
#' @return The upper bound on the boundary crossing probability for confidence sequences and corresponding K value.
# @export
# @examples
# cross_prob_cs(1, 1e+6)
# cross_prob_cs(2, 1e+8, 1e+2, eta = 2)
cross_prob_cs <- function(g, nmax,
                          nmin = 1L,
                          eta = 1.1) {
  if (nmin > nmax) stop("nmin must be smaller than or equal to nmax.")
  if (g <= 0) stop("g must be a positive number.")

  if (nmin == nmax){
    prob <- exp(-g)
    K <- 0L
    out <- list(prob = prob,
                K = K)
    return(out)
  }
  log_value <- log(nmax/nmin, base = eta)
  K <- ceiling(log_value)
  if (1 + log_value - K < 1e-14){
    K <- K - 1 # Correct possible floating point roundoff error.
  }

  if (nmin == 1){
    prob <- K * exp(-g / eta)
  } else {
    prob <- exp(-g) + K * exp(-g / eta)
  }
  out <- list(prob = prob,
              K = K)
  return(out)
}


#' Grid-searching the upper bound on the boundary crossing probability for confidence sequences.
#'
#' \code{cross_prob_cs_search} is used to grid-search the constant threshold for the GLR-like confidence sequences with finite target time interval.
#'
#' @param g Constant threshold (positive numeric).
#' @param nmax Upper bound of the target time interval.
#' @param nmin Lower bound of the target time interval. (default = 1L)
#' @param m_upper Upper bound on the grid for searching m value. (default = 1e+3L)

#' @return Grid-search result of the upper bound on the boundary crossing probability for confidence sequences.
# @export
# @examples
# cross_prob_cs_search(1, 1e+6)
# cross_prob_cs_search(2, 1e+8, 1e+2, m_upper = 500L)
cross_prob_cs_search <- function(g, nmax,
                                 nmin = 1L,
                                 m_upper = 1e+3L) {

  if (nmin > nmax) stop("nmin must be smaller than or equal to nmax.")
  if (g <= 0) stop("g must be a positive number.")

  if (nmin == nmax) {
    prob <- exp(-g)
    K <- 0L
    out <- list(prob = prob,
                K = K,
                eta = NA)
    return(out)
  }
  # m values for the grid-search
  m_vec <- seq(from = 1, to = m_upper)
  cal_eta <- function(m) (nmax / nmin)^(1/m)
  eta_vec <- sapply(m_vec, cal_eta)
  cal <- m_vec * exp(-g /eta_vec)

  if (nmin == 1){
    prob <- cal
  } else {
    prob <- exp(-g) + cal
  }
  out <- list(prob = prob,
              K = m_vec,
              eta = eta_vec)
  return(out)
}

#' Constant boundary value given a boundary crossing probability for confidence sequences.
#'
#' \code{const_boundary_cs} is used to compute the constant threshold for the GLR-like confidence sequences with finite target time interval.
#'
#' @param alpha An upper bound on the boundary crossing probability (positive numeric in \code{[1e-16,0.5]}).
#' @param nmax Upper bound of the target time interval. If both \code{nmax} and \code{d} are provided, \code{nmax} will be ignored. (default = NULL)
#' @param nmin Lower bound of the target time interval. (default = 1L)
#' @param d Bregman divergence between the null and alternative spaces. Either \code{d} or \code{nmax} must be specified (default = NULL).
#' @param m_upper Upper bound on the grid-search for m (default = 1e+3L).
#'
#' @return Constant threshold for GLR-like confidence sequence with level \code{alpha} based on Theorem 7.
#' @export
#' @examples
#' const_boundary_cs(0.05, 1e+6)
#' const_boundary_cs(0.025, 1e+8, 1e+2, m_upper = 100L)
const_boundary_cs <- function(alpha, nmax = NULL,
                              nmin = 1L,
                              d = NULL,
                              m_upper = 1e+3L) {
  if (alpha < 1e-16 | alpha > 0.5) stop("alpha must be in [1e-16,0.5].")

  if (is.null(nmax)){
    if (!is.null(d)){
      if (d <= 0) stop("d must be a positive number.")
      prob_diff <- function(g){
        nmax <- g / d
        min(cross_prob_cs_search(g, nmax, nmin, m_upper)$prob) - alpha
      }
    } else {
      stop("Either d (positive numeric) or nmax(positive integer) must be specified.")
    }
  } else {
    if (!is.null(d)) warning("Both d and nmax are provided - d will be ignored.")
    if (nmin > nmax) stop("nmin must be smaller than or equal to nmax.")
    prob_diff <- function(g){
      min(cross_prob_cs_search(g, nmax, nmin, m_upper)$prob) - alpha
    }
  }

  g <- stats::uniroot(prob_diff,
                      interval = c(log(1/alpha), 10 * log(1/alpha)),
                      extendInt = "downX",
                      tol = 1e-10)$root
  g_eta <- cross_prob_cs_search(g, nmax, nmin, m_upper)
  prob_min_ind <- which.min(g_eta$prob)
  prob <- g_eta$prob[prob_min_ind]
  eta <- g_eta$eta[prob_min_ind]
  K <- g_eta$K[prob_min_ind]

  out <- list(g = g,
              prob = prob,
              alpha = alpha,
              nmax = nmax,
              nmin = nmin,
              K = K,
              eta = eta
  )

  return(out)
}


