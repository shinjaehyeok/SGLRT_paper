test_that("const_boundary must yield accurate prob bound",{
  alpha <- 10^seq(-1,-15)
  n_vec <- 10^seq(1,6)
  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      nmax_vec <- n_vec[n_vec >= n]
      prob <- sapply(nmax_vec, function(n) const_boundary_cs(a, n, nmin)$prob)
      expect_true(max(abs(prob-a)) < 1e-12)
    }
  }
})

test_that("mixture_wrapper should yield the same CI as Gaussian mixture one",{

  const_gauss_mix_CI <- function(g, eta, K, nmax, nmin = 1L, sig = 1){
    # Stitching and discrete mixture for sub-Gaussian distributions
    # From constant boundary

    # Inputs psi_star_inv; d; eta; K; g
    psi_star_inv <- function(x){
      sig * sqrt(2 * x)
    }

    psi_star_derv <- function(x){
      x / sig^2
    }

    psi_star <- function(x){
      x^2 / 2 / sig^2
    }

    # Calculate weight (a / h)
    w <- exp(-g / eta)
    w_0 <- ifelse(nmin == 1, 0, exp(-g))

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
        w * sum(exp(v * (b_vec + a_vec * x_bar)))
      }
    } else if (nmin == nmax){
      dis_mart <- function(x_bar, v){
        w_0 * exp(v * (b_vec[1] + a_vec[1] * x_bar))
      }
    } else {
      dis_mart <- function(x_bar, v){
        w_0 * exp(v * (b_vec[1] + a_vec[1] * x_bar))  +
          w * sum(exp(v * (b_vec[-1] + a_vec[-1] * x_bar)))
      }
    }

    # Boundary function for discrete mixture
    dis_mix_fn <- function(v){
      f <- function(x) dis_mart(x, v = v) - 1
      upper <- stitch_fn(v)
      bound <- uniroot(f, c(upper / 2, upper * 1.001), tol = 1e-12)$root
      return(bound)
    }

    out <- list(stitch_fn = stitch_fn,
                dis_mix_fn = dis_mix_fn,
                sig = sig,
                nmax = nmax,
                nmin = nmin,
                g = g,
                eta = eta,
                K = K)
    return(out)
  }


  alpha <- 10^seq(-1,-15)
  n_vec <- 10^seq(1,6)
  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      nmax_vec <- n_vec[n_vec >= n]
      prob <- sapply(nmax_vec, function(n) const_boundary_cs(a, n, nmin)$prob)
      expect_true(max(abs(prob-a)) < 1e-12)
    }
  }
})

