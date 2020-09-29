# list of underlying bregman divergence or psi^* functions for each representative class

#' Pre-defined psi^* and Bregman divergence functions for sub-Gaussian family.
#'
#' @param sig The sigma parameter of the sub-Gaussian family (default = 1).
#' @param is_add If \code{is_add} is \code{TRUE} then return psi^* functions for \code{SGRL_CI_additive}. Otherwise, return Bregman divergence functions for \code{SGLR_CI}.
#'
#' @return A list of pre-defined psi^* and Bregman divergence functions for sub-Gaussian family.
#'
#' @export
#'
generate_sub_Gaussian_fn <- function(sig = 1, is_add = TRUE){
  if (is_add){
    G_add_fn_list <- list(psi_star_inv = function(x){sig * sqrt(2 * x)},
                          psi_star_derv = function(x){x / sig^2}
    )
    return(G_add_fn_list)
  } else {
    G_fn_list <- list(breg = function(mu_1, mu_0){(mu_1 - mu_0)^2 / 2 / sig^2},
                      breg_pos_inv = function(d, mu_0){mu_0 + sig * sqrt(2 * d)},
                      breg_neg_inv = function(d, mu_0){mu_0 - sig * sqrt(2 * d)},
                      breg_derv = function(z, mu_0){(z - mu_0) / sig^2},
                      mu_lower = NULL,
                      mu_upper = NULL,
                      grid_by = NULL)
    return(G_fn_list)
  }
}


#' Pre-defined Bregman divergence functions for sub-Bernulli family.
#'
#' @param grid_by The size of grid-window of mean space (default = 0.1).
#'
#' @return A list of pre-defined Bregman divergence functions for sub-Bernulli family.
#'
#' @export
generate_sub_ber_fn <- function(grid_by = 0.1){
  ber_fn_list <- list(breg = function(mu_1, mu_0){
    if (mu_1 == mu_0) return(0)
    if (mu_0 == 1 | mu_0 == 0){
      return(Inf)
    }
    if (mu_1 == 1){
      return(- log(mu_0))
    } else if (mu_1 == 0){
      return(- log(1 - mu_0))
    }
    mu_1 * log(mu_1 / mu_0) + (1-mu_1) * log((1-mu_1) / (1-mu_0))
  },
  breg_pos_inv = function(d, mu_0){
    if (mu_0 == 0) return(0)
    if (mu_0 == 1) return(1)
    f <- function(z) ber_fn_list$breg(z, mu_0) - d
    if (f(1) <= 0) return(1)
    z <- stats::uniroot(f,
                        c(mu_0, 1),
                        tol = 1e-12)$root
    return(z)
  },
  breg_neg_inv = function(d, mu_0){
    if (mu_0 == 0) return(0)
    if (mu_0 == 1) return(1)

    f <- function(z) ber_fn_list$breg(z, mu_0) - d
    if (f(0) <= 0) return(0)
    z <- stats::uniroot(f,
                        c(mu_0, 0),
                        tol = 1e-12)$root
    return(z)
  },
  breg_derv = function(z, mu_0){
    log(z * (1-mu_0) / (1-z) / mu_0)
  },
  mu_lower = 0,
  mu_upper = 1,
  grid_by = 0.1
  )
  return(ber_fn_list)
}
