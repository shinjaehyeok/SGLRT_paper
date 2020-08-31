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
