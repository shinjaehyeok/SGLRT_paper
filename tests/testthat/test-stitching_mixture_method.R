test_that("GLR-like function should yield a valid CS for sub Gaussian case",{
  alpha <- 10^seq(-1,-6)
  n_vec <- 10^seq(1,5)

  GLR_like_CI <- function(v, g, nmax, nmin = 1L){
    if (v < nmin){
      const_1 <- sqrt(2) + (nmin / v - 1) / sqrt(2)
      out <- sqrt(g / nmin) * const_1
    } else if (v > nmax){
      const_2 <- sqrt(2) + (nmax / v - 1) / sqrt(2)
      out <- sqrt(g / nmax) * const_2
    } else {
      out <- sqrt(2 * g / v)
    }
    return(out)
  }

  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)
      g_list <- const_boundary_cs(a, nmax, nmin)
      bound_GLR_direct <- sapply(n_vec, function(v) GLR_like_CI(v,
                                                                g_list$g,
                                                                nmax,
                                                                nmin))
      GLR_like_fn <- SGLR_CI(a, nmax, nmin)$GLR_like_fn
      bound_GLR_like <- sapply(n_vec, GLR_like_fn)

      expect_true(max(abs(bound_GLR_direct-bound_GLR_like)) < 1e-12)
    }
  }
})

test_that("GLR-like >= Stitching >= Discrete mixture",{
  alpha <- 10^seq(-1,-6)
  n_vec <- 10^seq(1,5)

  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)

      out <- SGLR_CI(a, nmax, nmin)
      v_vec <- sample(1e+6, 1e+2)
      bound_GLR_like <- sapply(v_vec, out$GLR_like_fn)
      bound_stitch <- sapply(v_vec, out$stitch_fn)
      bound_dis_mix <- sapply(v_vec, out$dis_mix_fn)

      expect_true(sum(bound_GLR_like < bound_stitch - 1e-12) == 0)
      expect_true(sum(bound_stitch < bound_dis_mix - 1e-12) == 0)
    }
  }
})

test_that("GLR-like >= Stitching >= Discrete mixture",{
  alpha <- c(0.1, 0.05, 0.01)
  n_vec <- 10^seq(1,4)
  B <- 1e+3
  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)

      out <- SGLR_CI(a, nmax, nmin)
      v_vec <- seq(1,1e+4)
      bound_GLR_like <- sapply(v_vec, out$GLR_like_fn)
      bound_stitch <- sapply(v_vec, out$stitch_fn)
      bound_dis_mix <- sapply(v_vec, out$dis_mix_fn)
      bound_chernoff <- sqrt(2 * log(1/a) / v_vec)

      errors <- data.frame(GLR = 0, stitch = 0, dis_mix = 0, point_chernoff = 0)

      for (i in 1:B){
        X_vec <- rnorm(1e+4)
        m_vec <- cumsum(X_vec) / v_vec
        errors$GLR <- errors$GLR +
          ifelse(sum(bound_GLR_like < m_vec)!=0, 1,0)
        errors$stitch <- errors$stitch +
          ifelse(sum(bound_stitch < m_vec)!=0, 1,0)
        errors$dis_mix <- errors$dis_mix +
          ifelse(sum(bound_dis_mix < m_vec)!=0, 1,0)
        errors$point_chernoff <- errors$point_chernoff +
          ifelse(sum(bound_chernoff < m_vec)!=0, 1,0)
      }
      error_rate <- errors / B
      expect_true(error_rate$GLR <= a + 2 * sqrt(a * (1-a) / B))
      expect_true(error_rate$stitch <= a + 2 * sqrt(a * (1-a) / B))
      expect_true(error_rate$dis_mix <= a + 2 * sqrt(a * (1-a) / B))
    }
  }
})
#
# plot(v_vec, m_vec, type = "l", ylim = c(-0.5, 0.5))
# lines(v_vec, m_vec + bound_GLR_like, col = 2)
# lines(v_vec, m_vec - bound_GLR_like, col = 2)
#
# lines(v_vec, m_vec + bound_stitch, col = 3)
# lines(v_vec, m_vec - bound_stitch, col = 3)
#
# lines(v_vec, m_vec + bound_dis_mix, col = 4)
# lines(v_vec, m_vec - bound_dis_mix, col = 4)
#
# lines(v_vec, m_vec + sqrt(2 * log(1/a) / v_vec), col = 5, lty = 2)
# lines(v_vec, m_vec - sqrt(2 * log(1/a) / v_vec), col = 5, lty = 2)
# abline(h = 0)
