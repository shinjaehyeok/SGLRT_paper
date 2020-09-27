#' GLR-like, Stitching and Discrete Mixture CS for additive sub-psi class
#'
#' \code{SGLR_CI_additive} is used to compute R functions which compute bounds of GLR-like, stitching and discrete mixture confidence sequences for additive sub-psi class designed to a finite target time interval.
#'
#' @inheritParams const_boundary_cs
#' @param psi_star_inv R function of \eqn{\psi_{+}^{*-1}} (Default: \eqn{\sqrt{2x}}) .
#' @param psi_star_derv R function of \eqn{\nabla\psi^*} (Default: \eqn{x}).
#' @param n_0 Lower bound of the sample size on which test statistics and CI will be computed (default = 1L).
#'
#' @return A list of R functions for GLR-like, stitching and discrete mixture bound which takes the sample size \code{n} as the input and return the distance from the sample mean to the upper bound of confidence interval at \code{n}. The list also contains related quantities to compute these bounds. See ADD_CITE for detailed explanations of these quantities.
#' \describe{
#'   \item{GLR_like_fn}{R function for the GLR-like bound.}
#'   \item{stitch_fn}{R function for the stitching bound.}
#'   \item{dis_mix_fn}{R function for the discrete mixture bound.}
#'   \item{log_mtg_fn}{R function for the log of the underlying log super-martingale.}
#'   \item{g}{The constant bouandry value.}
#'   \item{eta}{The eta value used to construct underlying martingales}
#'   \item{K}{The number of LR-like martingales used to construct bounds for \eqn{n \geq n_{\min}} part.}
#' }
#'
#' @export
#' @examples
#' SGLR_CI_additive(0.025, 1e+5)
#' SGLR_CI_additive(0.025, 1e+6, 1e+2)

SGLR_CI_additive <- function(alpha,
                             nmax,
                             nmin = 1L,
                             m_upper = 1e+3L,
                             psi_star_inv = function(x){sqrt(2 * x)},
                             psi_star_derv = function(x){x},
                             n_0 = 1L)
{
  # Calculate g, eta, K
  param_out <- const_boundary_cs(alpha, nmax, nmin, m_upper)
  g <- param_out$g
  eta <- param_out$eta
  K <- param_out$K

  # Compute GLR-like bounds
  psi_inv_val1 <- psi_star_inv(g / nmin)
  slop1 <- g / psi_star_derv(psi_inv_val1)
  const1 <- psi_inv_val1 - slop1 / nmin

  psi_inv_val2 <- psi_star_inv(g / nmax)
  slop2 <- g / psi_star_derv(psi_inv_val2)
  const2 <- psi_inv_val2 - slop2 / nmax


  GLR_like_fn <- function(v){
    if (v < n_0) return(Inf)
    if (v < nmin){
      out <- const1 + slop1 / v
    } else if (v > nmax){
      out <- const2 + slop2 / v
    } else {
      out <- psi_star_inv(g / v)
    }
    return(out)
  }


  # Compute stitching and discrete mixture bounds
  # Calculate a_vec (slop), b_vec (const)
  if (nmin == n_0) {
    g_eta_vec <-  g / eta^seq(1,K) / n_0
  } else {
    g_eta_vec <-  g / eta^seq(0,K) / nmin
  }
  inv_g_eta_vec <- sapply(g_eta_vec, psi_star_inv)
  a_vec <- sapply(inv_g_eta_vec, psi_star_derv)
  b_vec <- g_eta_vec - a_vec * inv_g_eta_vec

  # Calculate coeff (CI := (\bar{X} - \min_k (1/n s_k - c_k), \infty)
  if (nmin == n_0){
    s_vec <- g / eta / a_vec
  } else {
    s_vec <- c(g / a_vec[1], g / eta / a_vec[-1])
  }
  c_vec <- b_vec / a_vec

  # Boundary function for stitching
  stitch_fn <- function(v){
    if (v < n_0) return(Inf)
    out <- min(s_vec / v - c_vec)
    return(out)
  }

  # Discrete mixture
  if (nmin == nmax){
    dis_mart <- function(x_bar, v){
      if (v < n_0) return(-Inf)
      out <- v * (b_vec[1] + a_vec[1] * x_bar) - g
      return(out)
    }
  } else if (nmin == n_0){
    dis_mart <- function(x_bar, v){
      if (v < n_0) return(-Inf)
      inner <- v * (b_vec + a_vec * x_bar)
      inner_max <- max(inner)
      out <- inner_max + log(sum(exp(inner - inner_max))) - g / eta
      return(out)
    }

  } else {
    dis_mart <- function(x_bar, v){
      if (v < n_0) return(-Inf)
      inner <- v * (b_vec + a_vec * x_bar)
      inner[1] <- inner[1] - g
      inner[-1] <- inner[-1] - g / eta
      inner_max <- max(inner)
      out <- inner_max + log(sum(exp(inner - inner_max)))
    }
  }

  # Boundary function for discrete mixture
  dis_mix_fn <- function(v){
    if (v < n_0) return(Inf)
    f <- function(x) dis_mart(x, v = v)
    upper <- stitch_fn(v)
    bound <- stats::uniroot(f, c(upper / 2, upper * 1.001), tol = 1e-12)$root
    return(bound)
  }

  out <- list(GLR_like_fn = GLR_like_fn,
              stitch_fn = stitch_fn,
              dis_mix_fn = dis_mix_fn,
              log_mtg_fn = dis_mart,
              alpha = alpha,
              nmax = nmax,
              nmin = nmin,
              g = g,
              eta = eta,
              K = K,
              n_0 = n_0)
  return(out)
}


#' GLR-like and Discrete Mixture CS for general sub-psi class
#'
#' \code{SGLR_CI} is used to compute R functions which compute bounds of GLR-like, stitching and discrete mixture confidence sequences for general sub-psi class designed to a finite target time interval. For the additive sub-psi class, we recommend to use \code{SGLR_CI_additive} for more efficient computations.
#'
#' @inheritParams const_boundary_cs
#' @param breg R function of \eqn{D(\mu_1, \mu_0)} which takes \code{mu_1} and \code{mu_0} as inputs (Default: \eqn{(mu_1 - mu_0)^2 / 2}).
#' @param breg_pos_inv R function of inverse of the mapping \eqn{z \mapsto D(z, \mu_0):= d} on \eqn{z > \mu_0} which takes \code{d} and \code{mu_0} as inputs (Default: \eqn{\mu_0+\sqrt{2d}}).
#' @param breg_neg_inv R function of inverse of the mapping \eqn{z \mapsto D(z, \mu_0):= d} on \eqn{z < \mu_0} which takes \code{d} and \code{mu_0} as inputs (Default: \eqn{\mu_0 - \sqrt{2d}}).
#' @param breg_derv R function of \eqn{\nabla_z D(z, \mu_0)} which takes \code{z} and \code{mu_0} as inputs  (Default: \eqn{z - \mu_0}).
#' @param mu_lower Lower bound of the mean parameter space (default = NULL).
#' @param mu_upper Upper bound of the mean parameter space (default = NULL).
#' @param grid_by The size of grid-window of mean space. Default is \code{NULL}.
#' @param n_0 Lower bound of the sample size on which test statistics and CI will be computed (default = 1L).
#'
#' @return A list of R functions for GLR-like, stitching and discrete mixture bound which takes sample mean \code{x_bar} and sample size \code{n} as the input and return the anytime-valid confidence interval at \code{n}. The list also contains related quantities to compute these bounds. See ADD_CITE for detailed explanations of these quantities.
#' \describe{
#'   \item{GLR_like_fn}{R function for the GLR-like bound.}
#'   \item{dis_mix_fn}{R function for the discrete mixture bound.}
#'   \item{log_GLR_like_stat_generator}{R function to generate log of GLR-like statistic minus the threshold.}
#'   \item{log_dis_mart_generator}{R function to generate log of GLR-like statistic minus the threshold.}
#'   \item{alpha}{alpha valud used to construct GLR-like and discrete mixture bound functions}
#'   \item{g}{The boundary value for initial \code{nmin} and \code{nmax}.}
#'   \item{eta}{The eta value used to construct underlying martingales}
#'   \item{K}{The number of LR-like martingales used to construct bounds for \eqn{n \geq n_{\min}} part.}
#'   \item{mu_range}{Minimum and maximum of the mean parater space.}
#'   \item{grid_by}{The size of grid-window of mean space.}
#' }
#'
#' @export
#' @examples
#' SGLR_CI(0.025, 1e+5)
#' SGLR_CI(0.025, 1e+6, 1e+2)

SGLR_CI <- function(alpha,
                    nmax,
                    nmin = 1L,
                    m_upper = 1e+3L,
                    breg = function(mu_1, mu_0){(mu_1 - mu_0)^2 / 2},
                    breg_pos_inv = function(d, mu_0){mu_0 + sqrt(2 * d)},
                    breg_neg_inv = function(d, mu_0){mu_0 - sqrt(2 * d)},
                    breg_derv = function(z, mu_0){z - mu_0},
                    mu_lower = NULL,
                    mu_upper = NULL,
                    grid_by = NULL,
                    n_0 = 1L)
{
  # Calculate g, eta, K for common parameters for both generators
  param_out <- const_boundary_cs(alpha, nmax, nmin, m_upper)
  g <- param_out$g
  eta <- param_out$eta
  K <- param_out$K
  CI_grid <- NULL

  # Construct grid if grid size and mean range are provided.
  if (!is.null(grid_by)){
    if (is.null(mu_lower) | is.null(mu_upper)){
      warning("The range of mean space (mu_lower, mu_upper) cannot found. grid_by will be ignored.")
    } else {
      CI_grid <- seq(mu_lower, mu_upper, by = grid_by)
    }
  }
  # Initializing function for statistic generators
  initial_fn <- function(mu_0, is_pos){
    # Define  the trivial Statistic for bounded mean parameter case
    trivial_stat_fn <- function(x_bar, v){
      if (v < n_0) return(-Inf)
      out <- ifelse(x_bar == mu_0, -Inf, Inf)
      return(out)
    }

    # If mu_0 is at the boundary, the generator return a trivial statistics
    if (!is.null(mu_upper)){
      if (mu_0 == mu_upper){
        return(list(is_trivial = TRUE, stat_fn = trivial_stat_fn))
      }
    }
    if (!is.null(mu_lower)){
      if (mu_0 == mu_lower){
        return(list(is_trivial = TRUE, stat_fn = trivial_stat_fn))
      }
    }

    if (is_pos){ # Compute statistic for the right-sided case (mu_1 > mu_0)
      # Check whether nmin is large enough for the bounded mean space case.
      if (!is.null(mu_upper)){
        n_0 <- g / breg(mu_upper, mu_0)
        # If nmin is too small, update nmin and nmax.
        if (nmin <= n_0){
          nmin <- n_0
          if (nmax <= n_0){ # Due to the boundary case, we do not allow eff_min = nmin = nmax.
            nmin <- n_0 + 1
            nmax <- nmin
          }
          # Update boundary value for the updated nmin and nmax
          param_out <- const_boundary_cs(alpha, nmax, nmin, m_upper)
          g <- param_out$g
          eta <- param_out$eta
          K <- param_out$K
        }
      }
      # Define mean mapping for the left_sided case (mu_1 < mu_0)
      z_fn <- function(d) breg_pos_inv(d, mu_0)
    } else { # Compute statistic for the left_sided case (mu_1 < mu_0)
      # Check whether nmin is large enough for the bounded mean space case.
      if (!is.null(mu_lower)){
        n_0 <- g / breg(mu_lower, mu_0)
        if (nmin < n_0){
          nmin <- n_0
          nmax <- ifelse(nmax > nmin, nmax, nmin)
          param_out <- const_boundary_cs(alpha, nmax, nmin, m_upper)
          g <- param_out$g
          eta <- param_out$eta
          K <- param_out$K
        }
      }
      # Define mean mapping for the left_sided case (mu_1 < mu_0)
      z_fn <- function(d) breg_neg_inv(d, mu_0)
    }
    out <- list(is_trivial = FALSE,
                params = param_out,
                nmin = nmin,
                nmax = nmax,
                z_fn = z_fn,
                n_0_updated = n_0)
    return(out)
  }


  # Compute GLR-like statistics
  GLR_like_stat_generator <- function(mu_0, is_pos = TRUE){
    # Initialize parameters
    initial_out <- initial_fn(mu_0, is_pos)
    if (initial_out$is_trivial){
      out <- list(stat_fn = initial_out$stat_fn,
                  is_trivial = TRUE,
                  alpha = alpha,
                  nmax = nmax,
                  nmin = nmin,
                  g = g,
                  n_0 = n_0)
      return(out)
    }

    # Set parameters
    g <- initial_out$params$g
    eta <- initial_out$params$eta
    K <- initial_out$params$K
    nmin <- initial_out$nmin
    nmax <- initial_out$nmax
    n_0_updated <- initial_out$n_0_updated

    # Get the mean maping
    z_fn <- initial_out$z_fn

    # Compute GLR stat
    d2 <- g / nmin
    d1 <- g / nmax

    mu_2 <- z_fn(d2)
    slop_2 <- breg_derv(mu_2, mu_0)
    mu_1 <- z_fn(d1)
    slop_1 <- breg_derv(mu_1, mu_0)

    GLR_like_stat_fn <- function(x_bar, v){
      if (v <  n_0_updated) return(-Inf)
      if (v < nmin){
        out <- d2 + slop_2 * (x_bar - mu_2)
        out <- ifelse(is.nan(out), -Inf, out)
      } else if (v > nmax){
        mu_1 <- z_fn(d1)
        out <- d1 + slop_1 * (x_bar - mu_1)
        out <- ifelse(is.nan(out), -Inf, out)
      } else {
        if (is_pos & x_bar >= mu_0){
          out <- breg(x_bar, mu_0)
        } else if (!is_pos & x_bar <= mu_0){
          out <- breg(x_bar, mu_0)
        } else{
          out <- 0
        }
      }
      return(max(out,0) * v -g)
    }
    out <- list(stat_fn = GLR_like_stat_fn,
                is_trivial = FALSE,
                alpha = alpha,
                nmax = nmax,
                nmin = nmin,
                g = g,
                n_0 = n_0,
                n_0_updated = n_0_updated)
    return(out)
  }



  # Compute discrete mixture statistics
  dis_mart_generator <- function(mu_0, is_pos = TRUE){
    # Initialize parameters
    initial_out <- initial_fn(mu_0, is_pos)
    if (initial_out$is_trivial) {
      out <- list(stat_fn = initial_out$stat_fn,
                  is_trivial = TRUE,
                  alpha = alpha,
                  nmax = nmax,
                  nmin = nmin,
                  g = g,
                  n_0 = n_0)
      return(out)
    }

    # Set parameters
    g <- initial_out$params$g
    eta <- initial_out$params$eta
    K <- initial_out$params$K
    nmin <- initial_out$nmin
    nmax <- initial_out$nmax
    n_0_updated <- initial_out$n_0_updated

    # Get the mean maping
    z_fn <- initial_out$z_fn

    # Compute dis. mixture statistic
    if (nmin == nmax){ # Case 1. nmin = nmax: Use single LR-like stat
      d1 <- g / nmax
      z_val <- z_fn(d1)
      slop_val <- breg_derv(z_val, mu_0)

      dis_mart <- function(x_bar, v){
        if (v < n_0_updated) return(-Inf)
        out <- v * ( d1 + slop_val * (x_bar - z_val) ) - g
        out <- ifelse(is.nan(out), -Inf, out)
        return(out)
      }

    } else if (nmin == n_0) { # Case 2. 1 = nmin < nmax : Use K LR-like stat
      g_eta_vec <-  g / eta^seq(1,K)
      z_vec <- sapply(g_eta_vec, z_fn)
      slop_vec <- sapply(z_vec, function(z) breg_derv(z, mu_0))

      dis_mart <- function(x_bar, v){
        if (v < n_0_updated) return(-Inf)
        inner <- v * ( g_eta_vec + slop_vec * (x_bar - z_vec) )
        inner_max <- max(inner)
        out <- inner_max + log(sum(exp(inner - inner_max))) - g / eta
        out <- ifelse(is.nan(out), -Inf, out)
        return(out)
      }
    } else { # Case 3. 1 < nmin < nmax : add the base line LR-like to the  K LR-like,
      g_eta_vec <-  g / eta^seq(0,K) / nmin
      z_vec <- sapply(g_eta_vec, z_fn)
      slop_vec <- sapply(z_vec, function(z) breg_derv(z, mu_0))

      dis_mart <- function(x_bar, v){
        if (v < n_0_updated) return(-Inf)
        inner <- v * ( g_eta_vec + slop_vec * (x_bar - z_vec) )
        # inside_term <- ifelse(abs(x_bar - z_vec[1]) < 1e-12, g_eta_vec[1], slop_vec[1])
        # inner[1] <- v * inside_term - g
        inner[1] <- inner[1] - g
        inner[-1] <- inner[-1] - g / eta
        inner_max <- max(inner)
        out <- inner_max + log(sum(exp(inner - inner_max)))
        out <- ifelse(is.nan(out), -Inf, out)
        return(out)
      }
    }
    out <- list(stat_fn = dis_mart,
                is_trivial = FALSE,
                alpha = alpha,
                nmax = nmax,
                nmin = nmin,
                g = g,
                n_0 = n_0,
                n_0_updated = n_0_updated)
    return(out)
  }

  # Compute CI bound

  find_CI_fn <- function(stat_fn, x_bar, v, is_pos = TRUE){

    # Define the test function
    f_mu <- function(mu) stat_fn(x_bar, v, mu, is_pos)

    # Compute the multiply factor for the searching range.
    thres <- log(1 / alpha) / v
    if (v < nmin){
      m_factor <- max(100, 100 * log(nmin / v, base = 10))
    } else if (v > nmax){
      m_factor <- max(100, 100 * log(v / nmax, base = 10))
    } else{
      m_factor <- 5
    }

    # Define search range.
    if (is_pos){
      if (is.null(mu_lower)){
        # Use breg_inv functions for approximate searching range
        z_bound_with_sign <- breg_neg_inv(thres, x_bar) - x_bar
        search_range <- c(x_bar + m_factor * z_bound_with_sign,
                          x_bar + z_bound_with_sign / 3)
        if (f_mu(search_range[1]) <= 0) return(-Inf)
      } else {
        if (!is.null(CI_grid)){ # If mu_lower and grid are both provided, use grid-search to refine the search space.
          grid_n <- length(CI_grid)
          grid_inter <- CI_grid[CI_grid <= x_bar]
          stat_vec <- sapply(grid_inter, f_mu)
          ind <- which(stat_vec < 0)
          min_ind <- min(ind)
          if (min_ind == 1){
            return(grid_inter[1])
          } else {
            search_range <- grid_inter[c(min_ind - 1, min_ind)]
          }
        } else {# If mu_lower is provided but grid is not, use the naive search range
          search_range <- c(mu_lower, x_bar)
          if (f_mu(search_range[1]) <= 0) return(mu_lower)
        }
      }
    } else {
      if (is.null(mu_upper)){
        z_bound_with_sign <- breg_pos_inv(thres, x_bar) - x_bar
        search_range <- c(x_bar + z_bound_with_sign / 3,
                          x_bar + m_factor * z_bound_with_sign)
        if (f_mu(search_range[2]) <= 0) return(mu_upper)
      } else {
        if (!is.null(CI_grid)){ # If mu_lower and grid are both provided, use grid-search to refine the search space.
          grid_n <- length(CI_grid)
          grid_inter <- CI_grid[CI_grid >= x_bar]
          stat_vec <- sapply(grid_inter, f_mu)
          ind <- which(stat_vec < 0)
          max_ind <- max(ind)
          if (max_ind == length(grid_inter)){
            return(grid_inter[length(grid_inter)])
          } else {
            search_range <-grid_inter[c(max_ind, max_ind + 1)]
          }
        } else{ # If mu_upper is provided but grid is not, use the naive search range
          search_range <- c(x_bar, mu_upper)
          if (f_mu(search_range[2]) <= 0) return(mu_upper)
        }
      }
    }

    # Find the boundary of mean space
    bound <- stats::uniroot(f_mu,
                            search_range,
                            tol = 1e-12)$root

    return(bound)
  }

  # Compute GLR-like CI bound

  GLR_like_fn <- function(x_bar, v, is_pos = TRUE){
    if (v < n_0) {
      if (is_pos){
        out <- ifelse(is.null(mu_lower), -Inf, mu_lower)
      } else{
        out <- ifelse(is.null(mu_upper), Inf, mu_upper)
      }
      return(out)
    }
    GLR_like_inner <- function(x_bar, v, mu_0, is_pos){
      f <- GLR_like_stat_generator(mu_0, is_pos)$stat_fn
      return(f(x_bar, v))
    }
    find_CI_fn(GLR_like_inner , x_bar, v, is_pos)
  }

  # Compute dis. mixture CI bound
  dis_mix_fn <- function(x_bar, v, is_pos = TRUE){
    if (v < n_0) {
      if (is_pos){
        out <- ifelse(is.null(mu_lower), -Inf, mu_lower)
      } else{
        out <- ifelse(is.null(mu_upper), Inf, mu_upper)
      }
      return(out)
    }
    dis_mix_inner <- function(x_bar, v, mu_0, is_pos){
      f <- dis_mart_generator(mu_0, is_pos)$stat_fn
      return(f(x_bar, v))
    }
    find_CI_fn(dis_mix_inner, x_bar, v, is_pos)
  }

  out <- list(GLR_like_fn = GLR_like_fn,
              dis_mix_fn = dis_mix_fn,
              log_GLR_like_stat_generator = GLR_like_stat_generator,
              log_dis_mart_generator = dis_mart_generator,
              alpha = alpha,
              nmax = nmax,
              nmin = nmin,
              g = g,
              eta = eta,
              K = K,
              mu_range = c(mu_lower, mu_upper),
              grid_by = grid_by,
              n_0 = n_0)
  return(out)
}
