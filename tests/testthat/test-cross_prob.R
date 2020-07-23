test_that("Minimum of cross_prob_searchs over m spaces should be smaller than one from search over eta space", {
  # Testing g and d values
  alpha <- 10^seq(-1,-8)
  g_vec <- log(1/alpha)
  d_vec <- (alpha[1:4] * 100) ^2 / 2

  for (g in g_vec){
    for (d in d_vec){
      # cross_prob_search over m space
      prob1 <- min(cross_prob_search(g,d)$prob)

      # cross_prob_search over p_hat space
      prob2 <- min(cross_prob_search2(g,d)$prob)

      # print(c(prob1,prob2, prob2-prob1, g,d))

      expect_true(prob1 <= prob2)
    }
  }

  # Visual check
  # g <- g_vec[2]
  # d <- 0.5
  # out <-  cross_prob_search(g,d)
  # out2 <- cross_prob_search2(g,d)
  #
  # plot(out2$eta, log(out2$prob))
  # lines(out$eta, log(out$prob), type = "l")

})