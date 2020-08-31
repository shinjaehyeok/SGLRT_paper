#' Calculate the upper bound on the boundary crossing probability.
#'
#' Time-varying threshold for the non-separation case.
#'
#' @param g_fun positive function on \eqn{[1,\infty)} (R function object, default is \code{NULL}).
#' @param exp_g_fun positive function on \eqn{\mathbb{N} \times [1,\infty)} for the exponent term in the upper bound in Theorem 1 (R function object, default is \code{NULL}). Either \code{g_fun} or \code{exp_g_fun} must be non-null. If both are provided, \code{exp_g_fun} will be used.
#' @param eta Peeling size (default = 1.1).
#' @param k_max Maximum numbers of summations. (positive integer, default = 1e+4L)
#'
#' @return The upper bound on the boundary crossing probability based on Theorem 1.
#' @export
#' @examples
#' g_fun <- function(n) 5 + 2 * log(1 + log(n))
#' exp_g_fun <- function(k, eta) exp(-5 / eta) / (1 + k * log(eta))^(2/eta)
#' cross_prob_varying(g_fun, eta = 1.1)
#' cross_prob_varying(g_fun = NULL, exp_g_fun, eta = 2)
cross_prob_varying <- function(g_fun = NULL,
                               exp_g_fun = NULL,
                               eta = 1.1,
                               k_max = 1e+4L) {
  if (is.null(g_fun) & is.null(exp_g_fun)) {
    stop("Either g_fun or exp_g_fun must be a function.")
  } else if (!is.null(g_fun) & is.null(exp_g_fun)){
    f <- function(k) exp(-g_fun(eta^k) / eta)
  } else if (is.null(g_fun) & !is.null(exp_g_fun)){
    f <- function(k) exp_g_fun(k, eta)
  } else {
    warning("Both g_fun and exp_g_fun were provided. Only exp_g_fun is used.")
    f <- function(k) exp_g_fun(k, eta)
  }
  exp_terms <- sapply(1:k_max, f)
  if (exp_terms[k_max] > 1e-6) warning("k_max needs to be increased for more accurate calculation,")
  # prob = round(sum(exp_terms) + 0.001, digits = 3)
  out <- list(prob = sum(exp_terms),
              fun = f,
              eta = eta)
  return(out)
}


#' Calculate the upper bound on the boundary crossing probability by using a grid-search
#'
#' Time-varying threshold for the non-separation case.
#'
#' @param g_fun positive function on \eqn{[1,\infty)} (R function object, default is \code{NULL}).
#' @param exp_g_fun positive function on \eqn{\mathbb{N} \times [1,\infty)} for the exponent term in the upper bound in Theorem 1 (R function object, default is \code{NULL}). Either \code{g_fun} or \code{exp_g_fun} must be non-null. If both are provided, \code{exp_g_fun} will be used.
#' @param k_max Maximum numbers of summations. (positive integer, default = 1e+4L)
#' @param eta_upper Upper bound on the grid-search for eta (default = 2).
#' @param eta_grid_size Grid size for the grid-search for eta (default = 0.01).
#'
#' @return Grid-search result of the upper bound on the boundary crossing probability based on Theorem 1 for time-varying boundary.
#' @export
#' @examples
#' g_fun <- function(n) 5 + 2 * log(1 + log(n))
#' exp_g_fun <- function(k, eta) exp(-5 / eta) / (1 + k * log(eta))^(2/eta)
#' cross_prob_varying_search(g_fun, eta_grid_size = 0.1)
#' cross_prob_varying_search(g_fun = NULL, exp_g_fun)
cross_prob_varying_search <- function(g_fun = NULL,
                                      exp_g_fun = NULL,
                                      k_max = 1e+4L,
                                      eta_upper = 2,
                                      eta_grid_size = 0.01){
   if (eta_upper <= 1) stop("eta_upper must be larger than 1.")
  if (eta_upper < 1 + eta_grid_size) {
    warning("eta_grid_size is too large compared to eta_upper. eta_upper is used without line search")
  }

  # eta values for the grid-search
  eta_vec <- seq(from = 1 + eta_grid_size, eta_upper, by = eta_grid_size)


  f <- function(e) suppressWarnings(cross_prob_varying(g_fun, exp_g_fun, eta = e, k_max))

  cal <- sapply(eta_vec, f)
  out <- list(prob = unlist(cal["prob",]),
              eta = eta_vec)
  return(out)
}
