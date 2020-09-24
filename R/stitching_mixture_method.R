#' GLR-like, Stitching and Discrete Mixture CS for additive sub-psi class
#'
#' \code{GLR_stitch_mix_CI} is used to compute R functions which compute bounds of GLR-like, stitching and discrete mixture confidence sequences desgined to a finite target time interval.
#'
#' @inheritParams const_boundary_cs
#' @param psi_star R function of \eqn{\psi^*} (Default: \eqn{x^2 /2} (sub-Gussian with \eqn{\sigma = 1}.))
#' @param psi_star_inv R function of \eqn{\psi^{*-1}} (Default: \eqn{\sqrt{2x}}) .
#' @param psi_star_derv R function of \eqn{\nabla\psi^*} (Default: \eqn{x}).
#'
#' @return A list of R functions for stitching and discrete mixture bound which takes the sample size \code{n} as the input and return the distance from the sample mean to the upper bound of confidence interval at \code{n}. The list also contains related quantities to compute these bounds. See ADD_CITE for detailed explanations of these quantities.
#' \describe{
#'   \item{GLR_like_fn}{R function for the GLR-like bound.}
#'   \item{stitch_fn}{R function for the stitching bound.}
#'   \item{dis_mix_fn}{R function for the discrete mixture bound.}
#'   \item{g}{The constant bouandry value.}
#'   \item{eta}{The eta value used to construct underlying martingales}
#'   \item{K}{The number of LR-like martingales used to construct bounds for \eqn{n \geq n_{\min}} part.}
#' }
#'
#' @export
#' @examples
#' GLR_stitch_mix_CI(0.025, 1e+5)
#' GLR_stitch_mix_CI(0.025, 1e+6, 1e+2)

GLR_stitch_mix_CI <- function(alpha,
                          nmax,
                          nmin = 1L,
                          m_upper = 1e+3L,
                          psi_star = function(x){x^2 / 2},
                          psi_star_inv = function(x){sqrt(2 * x)},
                          psi_star_derv = function(x){x})
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
  if (nmin == 1) {
    g_eta_vec <-  g / eta^seq(1,K)
  } else {
    g_eta_vec <-  g / eta^seq(0,K) / nmin
  }
  inv_g_eta_vec <- sapply(g_eta_vec, psi_star_inv)
  a_vec <- sapply(inv_g_eta_vec, psi_star_derv)
  b_vec <- g_eta_vec - a_vec * inv_g_eta_vec

  # Calculate coeff (CI := (\bar{X} - \min_k (1/n s_k - c_k), \infty)
  if (nmin == 1){
    s_vec <- g / eta / a_vec
  } else {
    s_vec <- c(g / a_vec[1], g / eta / a_vec[-1])
  }
  c_vec <- b_vec / a_vec

  # Boundary function for stitching
  stitch_fn <- function(v){
    min(s_vec / v - c_vec)
  }

  # Discrete mixture
  if (nmin == 1){
    dis_mart <- function(x_bar, v){
      inner <- v * (b_vec + a_vec * x_bar)
      inner_max <- max(inner)
      out <- inner_max + log(sum(exp(inner - inner_max))) - g / eta
      return(out)
    }
  } else if (nmin == nmax){
    dis_mart <- function(x_bar, v){
      v * (b_vec[1] + a_vec[1] * x_bar) - g
    }
  } else {
    dis_mart <- function(x_bar, v){
      inner <- v * (b_vec + a_vec * x_bar)
      inner[1] <- inner[1] - g
      inner[-1] <- inner[-1] - g / eta
      inner_max <- max(inner)
      out <- inner_max + log(sum(exp(inner - inner_max)))
    }
  }

  # Boundary function for discrete mixture
  dis_mix_fn <- function(v){
    f <- function(x) dis_mart(x, v = v)
    upper <- stitch_fn(v)
    bound <- stats::uniroot(f, c(upper / 2, upper * 1.001), tol = 1e-12)$root
    return(bound)
  }

  out <- list(GLR_like_fn = GLR_like_fn,
              stitch_fn = stitch_fn,
              dis_mix_fn = dis_mix_fn,
              alpha = alpha,
              nmax = nmax,
              nmin = nmin,
              g = g,
              eta = eta,
              K = K)
  return(out)
}

