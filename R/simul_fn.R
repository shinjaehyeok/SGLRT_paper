# Conduct a simulation to compare sequential tests and fixed sample size test.
#'
#' \code{test_simul} is used to conduct a simulation which compare the sample size and detection rates of sequential GLR-lik, discrete mixture test and a fixed sample size test based on the sample mean.
#'
#' @param mu_target The boundary of the target alternative hypothesis space.
#' @param mu_true The underlying mean of the data-generating distribution (default = \code{mu_target}).
#' @param mu_0 The boundary of the null hypothesis space.
#' @param alpha The upper bound of the type-I error (default = \code{0.1})
#' @param beta The upper bound of the type-II error (default = \code{0.1}).
#' @param B The number of the repeats of the simulation (default = \code{100}).
#' @param print_progress A logical value. If \code{print_progress} = \code{TRUE}, a progress bar of the simulation will be printed (default = \code{TRUE}).
#' @param print_result A logical value. If \code{print_result} = \code{TRUE}, A short summary of the simulation result will be presented (default = \code{FALSE}).
#' @param sample_generator R function to generate i.i.d. random samples. The function must take the number of samples \code{n} and the underlying mean \code{mu_true} as inputs.
#' @param power_fn R function to compute the power of the test given the sample size. The function must take the number of samples \code{n}, the boundary of the alternative space \code{mu_target}, the boundary of the null space \code{mu_0}, and type-I error bound \code{alpha}.
#' @param fixed_test_fn R function to conduct a fixed sample size test. The function must take the sample mean \code{x_bar}, the number of sample \code{n}, The boundary of the null space \code{mu_0}, and the type-I error bound \code{alpha}. The function must return a positive numeric value if the test reject the null.
#' @param seq_test_fn R function to conduct sequential GLR-like and discrete mixture tests. The function must take the type-I error bound \code{alpha}, the upper bound \code{nmax} and the lower bound \code{nmin} of the target interval.
#'
#' @return A list of simulation results and underlying settings.
#' \describe{
#'   \item{reject_rate}{Estimated probabilities of rejecting the null hypothesis.}
#'   \item{sample_size}{Estimated average sample sizes of testing procedures.}
#'   \item{early_stop_ratio}{Estimated probabilities of tests being stopped earlier than the fixed sample size test.}
#'   \item{reject_list}{A list of vectors representing whether each test reject the null (\code{1}) or not (\code{0}) for each simulation run.}
#'   \item{sample_size_list}{A list of numbers of samples used in each test.}
#'   \item{mu_target}{The boundary of the target alternative hypothesis space.}
#'   \item{mu_true}{The underlying mean of the data-generating distribution.}
#'   \item{mu_0}{The boundary of the null hypothesis space.}
#'   \item{alpha}{The upper bound of the type-I error.}
#'   \item{beta}{The upper bound of the type-II error.}
#' }
#'
#' @export
#' @examples
#' # For the example of the simulation, please check https://github.com/shinjaehyeok/SGLRT_paper.
test_simul <- function(mu_target,
                       mu_true = mu_target,
                       mu_0,
                       alpha = 0.1,
                       beta = 0.1,
                       B = 100,
                       print_progress = TRUE,
                       print_result = FALSE,
                       sample_generator = G_sample,
                       power_fn = G_pwr,
                       fixed_test_fn = z_test,
                       seq_test_fn = seq_G_test_generator){

  # Calculate sample size to achieve the target power by the fixed sample size test
  if (mu_target != mu_0){
    sample_size <- stats::uniroot(function(n){
      power_fn(n, mu_target, mu_0, alpha) - 1 + beta
    },
    c(10, 1e+4))$root
    n_exact <- ceiling(sample_size)
  } else {
    n_exact <- 1e+3
  }

  is_pos <- ifelse(mu_target >= mu_0, 1, 0)

  # Set terminal time of sequential test.
  n <- 2 * n_exact

  # Set target time interval
  nmax <- n
  nmin <- round(max(10, n_exact / 10))

  # Range of experimentation
  n_vec <- seq(1, n)

  # Define test functions
  test_list <- list(exact = NULL,
                    GLR_like = NULL,
                    dis_mix = NULL
  )

  # Exact test
  test_list$exact <- function(x_bar, n){
    exact_out <-  fixed_test_fn(x_bar, n, mu_0, alpha)
    return(exact_out)
  }


  # Sequential tests
  G_out <- seq_test_fn(alpha, nmax, nmin)


  GLR_test <- G_out$log_GLR_like_stat_generator(mu_0, is_pos)
  test_list$GLR_like <- GLR_test$stat_fn

  dis_mix_test <- G_out$log_dis_mart_generator(mu_0, is_pos)
  test_list$dis_mix <- dis_mix_test$stat_fn


  # Type 1 error or power simulation
  is_reject <- function(test_value){
    num_reject <- sum(test_value > 0)
    return(ifelse(num_reject > 0, 1 ,0))
  }

  reject <- list(exact_hacking = numeric(B),
                 GLR_like = numeric(B),
                 dis_mix = numeric(B),
                 exact_test = numeric(B))

  sample_size_list <- list(exact_hacking = numeric(B),
                           GLR_like = numeric(B),
                           dis_mix = numeric(B),
                           exact_test = numeric(B))
  i <- 0
  if (print_progress)  pb <- utils::txtProgressBar(min = 0, max = B, style = 3)
  while (i < B){
    i <- i + 1
    x_vec <- sample_generator(n, mu_true)
    x_bar_vec <- cumsum(x_vec) / n_vec
    test_value_list <- list(exact_hacking = NULL,
                            GLR_like = NULL,
                            dis_mix = NULL,
                            exact_test = NULL
    )
    for (j in seq_along(test_list)){
      test_value_list[[j]] <- mapply(test_list[[j]], x_bar_vec, n_vec)
    }
    test_value_list$exact_test <-  test_value_list$exact_hacking[n_exact]

    for (k in seq_along(test_value_list)){
      if (k == length(test_value_list)) next
      reject_n <- which(test_value_list[[k]] >= 0)
      stopped_n <- ifelse(length(reject_n) == 0, n, min(reject_n))
      sample_size_list[[k]][i] <- stopped_n
      reject[[k]][i] <- is_reject(test_value_list[[k]])
    }
    sample_size_list$exact_test[i] <- n_exact
    reject$exact_test[i] <- is_reject(test_value_list$exact_test)
    if (print_progress) utils::setTxtProgressBar(pb, i)
  }
  if (print_progress)  close(pb)

  reject_rate <- sapply(reject, mean)
  sample_size <- sapply(sample_size_list, mean)
  early_stop_ratio <- sapply(seq(1, length(sample_size_list)),
                             function(k) mean(sample_size_list[[k]] <= n_exact))
  out <- list(reject_rate = reject_rate,
              sample_size = sample_size,
              early_stop_ratio = early_stop_ratio,
              reject_list = reject,
              sample_size_list = sample_size_list,
              mu_target = mu_target,
              mu_true = mu_true,
              mu_0 = mu_0,
              alpha = alpha,
              beta = beta)

  if (print_result){
    print(reject_rate)
    print(sample_size)
    print(early_stop_ratio)
    graphics::hist(sample_size_list$dis_mix, xlim = c(1, n))
    graphics::abline(v = n_exact)
  }
  return(out)
}


# Functions for tests based on a fixed sample size.

# Compute power of Z-test
G_pwr <- function(n, mu_target, mu_0, alpha){
  thres_exact <- stats::qnorm(alpha, mu_0, sd = 1/sqrt(n), lower.tail = FALSE)
  pwr <- stats::pnorm(thres_exact, mu_target, sd = 1/sqrt(n), lower.tail = FALSE)
  return(pwr)
}

# Z-test function (reject the null if the returned value is positive)
z_test <- function(x_bar, n, mu_0, alpha){
  z_alpha <- stats::qnorm(1-alpha)
  x_bar - mu_0 - z_alpha / sqrt(n)
}

# Compute power of binomial test
binom_pwr <- function(n, mu_target, mu_0, alpha){
  n <- floor(n)
  thres_exact <- stats::qbinom(alpha, n, mu_0, lower.tail = FALSE)
  pwr <- stats::pbinom(thres_exact, n, mu_target, lower.tail = FALSE)
  return(pwr)
}

# Exact binomial test (reject the null if the returned value is positive)
binom_test <- function(x_bar, n, mu_0, alpha){
  thres_exact <- stats::qbinom(alpha, n, mu_0, lower.tail = FALSE)
  return(x_bar * n - thres_exact)
}

# Sequential tests

# GLR-like and discrete mixture test for sub-Gaussian
seq_G_test_generator <- function(alpha, nmax, nmin){
  G_out <- SGLR_CI(alpha, nmax, nmin)
  return(G_out)
}

# GLR-like and discrete mixture test for sub-Bernoulli
ber_fn_list <- generate_sub_ber_fn()
seq_ber_test_generator <- function(alpha, nmax, nmin){
  ber_out <- SGLR_CI(alpha,
                     nmax,
                     nmin,
                     breg = ber_fn_list$breg,
                     breg_pos_inv = ber_fn_list$breg_pos_inv,
                     breg_neg_inv = ber_fn_list$breg_neg_inv,
                     breg_derv = ber_fn_list$breg_derv,
                     mu_lower = ber_fn_list$mu_lower,
                     mu_upper = ber_fn_list$mu_upper,
                     grid_by = ber_fn_list$grid_by)
  return(ber_out)
}

# Sample generators
# Gaussian samples
G_sample <- function(n, mu_true){
  stats::rnorm(n, mean = mu_true)
}

# Bernoulli samples
ber_sample <- function(n, mu_true){
  stats::rbinom(n, 1, prob = mu_true)
}

