#' Calculate the upper bound on the boundary crossing probability.
#'
#' Constant threshold for well-seperated case.
#'
#' @param g Constant threshold (positive numeric).
#' @param d Bregman divergence between the null and alternative spaces (positive numeric).
#' @param eta Peeling size (default = 1.36).
#'
#' @return The upper bound on the boundary crossing probability based on Theorem 1 and corresponding K value.
#' @export
#' @examples
#' cross_prob(1, 1)
#' cross_prob(2, 1, eta = 2)
cross_prob <- function(g, d,
                       eta = 1.36) {
  if (min(g, d) <= 0) stop("g and d must be a positive number.")
  if (g/d <= 1) {
    prob <- exp(-g)
    K <- 1
  } else {
    K <- ceiling(log(g/d, base = eta))
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
#' @param m_upper Upper bound on the grid-search for m (default = 1e+3).

#' @return Grid-search result of the upper bound on the boundary crossing probability based on Theorem 1.
#' @export
#' @examples
#' cross_prob_search(1, 1)
#' cross_prob_search(2, 1, m_upper = 500)
cross_prob_search <- function(g, d,
                       m_upper = 1e+3) {

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
