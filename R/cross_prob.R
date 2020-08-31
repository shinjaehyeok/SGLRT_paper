#' Calculate the upper bound on the boundary crossing probability.
#'
#' Constant threshold for well-seperated case.
#'
#' @param g Constant threshold (positive numeric).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#' @param eta Peeling size (default = 1.1).
#'
#' @return The upper bound on the boundary crossing probability based on Theorem 1 and corresponding K value.
#' @export
#' @examples
#' cross_prob(1, 1)
#' cross_prob(2, 1, eta = 2)
cross_prob <- function(g, d,
                       eta = 1.1) {
  if (min(g, d) <= 0) stop("g and d must be a positive number.")
  if (g/d <= 1) {
    prob <- exp(-g)
    K <- 1
  } else {
    log_value <- log(g/d, base = eta)
    K <- ceiling(log_value)
    if (1 + log_value - K < 1e-14){
      K <- K - 1 # Correct possible floating point roundoff error.
    }
    prob <- K * exp(-g / eta)
  }
  out <- list(prob = prob,
              K = K)
  return(out)
}

#' Calculate the upper bound on the boundary crossing probability by using      Lorden's formula.
#'
#' Constant threshold for well-seperated case.
#'
#' @inheritParams cross_prob
#'
#' @return The upper bound on the boundary crossing probability based on Lorden's formula
#' @export
#'
#' @examples
#' cross_prob_lorden(1, 1)
#' cross_prob_lorden(1, 2)
cross_prob_lorden <- function(g, d){
  if (g / d <= 1) {
    out <- exp(-g)
  } else {
    out <- (1 + g/d) * exp(-g)
  }
  return(out)
}


#' Calculate the upper bound on the boundary crossing probability by using a grid-search
#'
#' Constant threshold for well-seperated case.
#'
#' @param g Constant threshold (positive numeric).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#' @param m_upper Upper bound on the grid-search for m (default = 1e+3L).

#' @return Grid-search result of the upper bound on the boundary crossing probability based on Theorem 1.
#' @export
#' @examples
#' cross_prob_search(1, 1)
#' cross_prob_search(2, 1, m_upper = 500L)
cross_prob_search <- function(g, d,
                       m_upper = 1e+3L) {

  if (g/d <= 1) {
    prob <- exp(-g)
    K <- 1
    out <- list(prob = prob,
                K = K,
                eta = NA)
    return(out)
  }
  # m values for the grid-search
  m_vec <- seq(from = 1, to = m_upper)
  cal_eta <- function(m) (g/d)^(1/m)
  eta_vec <- sapply(m_vec, cal_eta)
  cal <- m_vec * exp(-g /eta_vec)
  # f <- function(m) m * exp(-g * (d/g)^(1/m))
  # cal <- sapply(m_vec, f)
  out <- list(prob = cal,
              K = m_vec,
              eta = eta_vec)
  return(out)
}

#' Calculate the constant boundary value given a boundary crossing probability.
#'
#' Constant threshold for well-seperated case.
#'
#' @param alpha An upper bound on the boundary crossing probability (positive numeric in \code{[1e-16,0.5]}).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#' @param m_upper Upper bound on the grid-search for m (default = 1e+3L).
#'
#' @return Constant threshold for GLR-lik statistic which makes the boundary crossing probability upper bounded by \code{alpha} based on Theorem 1.
#' @export
#' @examples
#' const_boundary(0.05, .1)
#' const_boundary(0.025, 1e-3, m_upper = 100L)
const_boundary <- function(alpha, d,
                           m_upper = 1e+3L) {
  if (alpha < 1e-16 | alpha > 0.5) stop("alpha must be in [1e-16,0.5].")
  if (d <= 0) stop("d must be a positive number.")
  prob_diff <- function(g){
    min(cross_prob_search(g, d, m_upper)$prob) - alpha
  }
  g <- stats::uniroot(prob_diff,
          interval = c(log(1/alpha), 10 * log(1/alpha)),
          extendInt = "downX",
          tol = 1e-10)$root
  g_eta <- cross_prob_search(g, d, m_upper)
  prob_min_ind <- which.min(g_eta$prob)
  prob <- g_eta$prob[prob_min_ind]
  eta <- g_eta$eta[prob_min_ind]
  K <- g_eta$K[prob_min_ind]

  out <- list(g = g,
              prob = prob,
              alpha = alpha,
              d = d,
              K = K,
              eta = eta
             )

  return(out)
}


#' Calculate the constant boundary value given a boundary crossing probability by using Lorden's bound.
#'
#' Constant threshold for well-seperated case.
#'
#' @param alpha An upper bound on the boundary crossing probability (positive numeric in \code{[1e-16,0.5]}).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#'
#' @return Constant threshold for GLR-lik statistic which makes the boundary crossing probability upper bounded by \code{alpha} based on the Lorden's bound.
#' @export
#' @examples
#' const_boundary_lorden(0.05, .1)
#' const_boundary_lorden(0.025, 1e-3)
const_boundary_lorden <- function(alpha, d) {
  if (alpha < 1e-16 | alpha > 0.5) stop("alpha must be in [1e-16,0.5].")
  if (d <= 0) stop("d must be a positive number.")
  prob_diff <- function(g){
    cross_prob_lorden(g,d) - alpha
  }
  g <- stats::uniroot(prob_diff,
               interval = c(log(1/alpha), 20 * log(1/alpha)),
               extendInt = "downX",
               tol = 1e-10)$root
  prob <- cross_prob_lorden(g, d)
  out <- list(g = g,
              prob = prob,
              alpha = alpha,
              d = d
              )

  return(out)
}


#' Calculate the upper bound on the boundary crossing probability by using a grid-search
#'
#' Constant threshold for well-seperated case.
#'
#' @param g Constant threshold (positive numeric).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#' @param eta_upper Upper bound on the grid-search for eta (default = g/d).
#' @param eta_grid_size Grid size for the grid-search for eta (default = 0.05).

#' @return Grid-search result of the upper bound on the boundary crossing probability based on Theorem 1.
#' @export
#' @examples
#' cross_prob_search2(1, 1)
#' cross_prob_search2(2, 1, eta_upper = 5, eta_grid_size = 0.5)
cross_prob_search2 <- function(g, d,
                              eta_upper = g/d,
                              eta_grid_size = 0.05) {
  if (g/d <= 1) {
    out <- list(prob = exp(-g),
                K = 1,
                eta = NA)
    return(out)
  }
  if (eta_upper <= 1) stop("eta_upper must be larger than 1.")
  if (eta_upper < 1 + eta_grid_size) warning("eta_grid_size is too large compared to eta_upper. eta_upper is used without line search")
  if (eta_upper >= 1e+3) eta_upper <- 1e+3

  # eta values for the grid-search
  eta_vec <- seq(from = 1 + eta_grid_size, eta_upper, by = eta_grid_size)


  f <- function(e) cross_prob(g, d, eta = e)
  cal <- sapply(eta_vec, f)
  out <- list(prob = unlist(cal["prob",]),
              K = unlist(cal["K",]),
              eta = eta_vec)
  return(out)
}
