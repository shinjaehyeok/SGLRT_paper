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
      GLR_like_fn <- SGLR_CI_additive  (a, nmax, nmin)$GLR_like_fn
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

      out <- SGLR_CI_additive  (a, nmax, nmin)
      v_vec <- sample(1e+6, 1e+2)
      bound_GLR_like <- sapply(v_vec, out$GLR_like_fn)
      bound_stitch <- sapply(v_vec, out$stitch_fn)
      bound_dis_mix <- sapply(v_vec, out$dis_mix_fn)

      expect_true(sum(bound_GLR_like < bound_stitch - 1e-12) == 0)
      expect_true(sum(bound_stitch < bound_dis_mix - 1e-12) == 0)
    }
  }
})

test_that("type1 error control",{
  alpha <- c(0.1, 0.05, 0.01)
  n_vec <- 10^seq(1,4)
  B <- 1e+3
  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)

      out <- SGLR_CI_additive(a, nmax, nmin)
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


test_that("SGLR_CI and SGLR_CI_additive must provide same bounds for sub Gaussian",{
  alpha <- 10^seq(-1,-6)
  n_vec <- 10^seq(1,5)

  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)

      out_add <- SGLR_CI_additive(a, nmax, nmin)
      out <- SGLR_CI(a, nmax, nmin)

      v_vec <- 10^(seq(0,5 + 0.2, by = 0.1))
      bound_GLR_like_add <- sapply(v_vec, out_add$GLR_like_fn) * sqrt(v_vec)
      bound_GLR_like <- sapply(v_vec, function(v) out$GLR_like_fn(0,v,is_pos = F))* sqrt(v_vec)
      bound_GLR_like2 <- -sapply(v_vec, function(v) out$GLR_like_fn(0,v,is_pos = T))* sqrt(v_vec)


      expect_true(max(abs(bound_GLR_like - bound_GLR_like_add)) < 1e-8)
      expect_true(max(abs(bound_GLR_like2 - bound_GLR_like_add)) < 1e-8)

      bound_dis_mix_add <- sapply(v_vec, out_add$dis_mix_fn) * sqrt(v_vec)
      bound_dis_mix <- sapply(v_vec, function(v) out$dis_mix_fn(0,v,is_pos = F))* sqrt(v_vec)
      bound_dis_mix2 <- -sapply(v_vec, function(v) out$dis_mix_fn(0,v,is_pos = T))* sqrt(v_vec)

      expect_true(max(abs(bound_dis_mix - bound_dis_mix_add)) < 1e-8)
      expect_true(max(abs(bound_dis_mix2 - bound_dis_mix_add)) < 1e-8)

    }
  }
})

test_that("Grid method and binary search must be consistent",{
  alpha <- 10^seq(-1,-3)
  n_vec <- 10^seq(1,5)

  for (a in alpha){
    for (n in n_vec){
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)

      out_grid <- SGLR_CI(a, nmax, nmin,CI_grid = seq(-100,100, by = 1))
      out <- SGLR_CI(a, nmax, nmin)

      v_vec <- 10^(seq(1,5 + 0.2, by = 0.1))
      bound_GLR_like <- sapply(v_vec,
                               function(v) out$GLR_like_fn(0,v,is_pos = F))
      bound_GLR_like_grid <- sapply(v_vec,
                                    function(v) out_grid$GLR_like_fn(0,v,is_pos = F))
      expect_true(sum(bound_GLR_like > bound_GLR_like_grid) == 0 )
      expect_true(sum(bound_GLR_like + 1 < bound_GLR_like_grid) == 0 )

      bound_dis_mix <- sapply(v_vec,
                              function(v) out$dis_mix_fn(0,v,is_pos = F))
      bound_dis_mix_grid <- sapply(v_vec,
                                   function(v) out_grid$dis_mix_fn(0,v,is_pos = F))
      expect_true(sum(bound_dis_mix > bound_dis_mix_grid) == 0 )
      expect_true(sum(bound_dis_mix + 1 < bound_dis_mix_grid) == 0 )
      # print(c(a, nmin, nmax))
    }
  }
})

test_that("Grid method and binary search must be consistent for sub-Bernouill",{

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
  })

  alpha <- 10^seq(-1,-3)
  n_vec <- 10^seq(1,5)


  for (a in alpha){
    for (n in n_vec){
      x_bar <- sample(c(0,0.1,0.5,0.9,1), 1)
      nmin <- n
      large_ind <- which(n_vec >= n)
      nmax <- ifelse(length(large_ind) > 1, sample(n_vec[n_vec >= n], 1), n)
      out_grid <- SGLR_CI(a, nmax, nmin,
                          breg = ber_fn_list$breg,
                          breg_pos_inv = ber_fn_list$breg_pos_inv,
                          breg_neg_inv = ber_fn_list$breg_neg_inv,
                          breg_derv = ber_fn_list$breg_derv,
                          CI_grid = seq(0,1, by =0.01),
                          mu_lower = 0,
                          mu_upper = 1)

      out <- SGLR_CI(a, nmax, nmin,
                     breg = ber_fn_list$breg,
                     breg_pos_inv = ber_fn_list$breg_pos_inv,
                     breg_neg_inv = ber_fn_list$breg_neg_inv,
                     breg_derv = ber_fn_list$breg_derv,
                     CI_grid = NULL,
                     mu_lower = 0,
                     mu_upper = 1)

      v_vec <- c(seq(1,9), round(10^(seq(1,5 + 0.2, by = 0.1))))

      bound_GLR_like <- sapply(v_vec,
                               function(v) out$GLR_like_fn(x_bar,v,is_pos = F))
      bound_GLR_like_grid <- sapply(v_vec,
                                    function(v) out_grid$GLR_like_fn(x_bar,v,is_pos = F))

      # plot(v_vec, bound_GLR_like_grid, log = "x", type = "l")
      # lines(v_vec, bound_GLR_like, col = 2)

      expect_true(sum(bound_GLR_like > bound_GLR_like_grid) == 0 )
      expect_true(sum(bound_GLR_like + .01 < bound_GLR_like_grid) == 0 )

      bound_dis_mix <- sapply(v_vec,
                               function(v) out$dis_mix_fn(x_bar,v,is_pos = F))
      bound_dis_mix_grid <- sapply(v_vec,
                                    function(v) out_grid$dis_mix_fn(x_bar,v,is_pos = F))




      # plot(v_vec, bound_dis_mix_grid, log = "x")
      # points(v_vec, bound_dis_mix, col = 2)

      expect_true(sum(bound_dis_mix > bound_dis_mix_grid) == 0 )
      expect_true(sum(bound_dis_mix + .01 < bound_dis_mix_grid) == 0 )

      # plot(v_vec, bound_GLR_like, log = "x")
      # points(v_vec, bound_dis_mix, col = 2)
      expect_true(sum(bound_dis_mix > bound_GLR_like + 1e-8) == 0 )


      plot(v_vec, bound_GLR_like_grid, log = "x", type = "l",
           main = paste0( c(nmin, nmax, x_bar, a)))
      lines(v_vec, bound_GLR_like, col = 2)
      lines(v_vec, bound_dis_mix_grid, col = 3)
      lines(v_vec, bound_dis_mix, col = 4)
    }
  }
})

# "sub-Bernouill must yield narrower CI than sub-G for bounded RVs"
