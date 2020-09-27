
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SGLRT

<!-- badges: start -->

<!-- badges: end -->

`SGLRT` is a R package implementation of Sequential Generalized
Likelihood Ratio (GLR)-like Tests and confidence seqeunces in **Paper**.

## Installation

You can install the `SGLRT` package from [GitHub](https://github.com/)
with:

``` r
# install.packages("devtools")
devtools::install_github("shinjaehyeok/SGLRT_public")
```

<!---You can install the released version of SGLRT from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("SGLRT")
```
-->

## How to reproduce all plots and simulation results in **paper**.

To reproduce plots, you need to install `latex2exp` package which parses
and converts LaTeX math formulas to R’s plotmath expressions. It is is
not installed, you can run the following command to install is.

``` r
install.packages("latex2exp")
```

### 1\. Compare boundaries for sequential GLR tests. *(Fig. 3 in Section III-A)*

For a given exponential family distribution, let
\(\mathrm{GLR}_n (> \mu_1, \leq\mu_0)\) be the GLR statistic based \(n\)
i.i.d. observations for the one-sided testing problem: \[
H_0 : \mu \leq \mu_0~~\text{vs}~~H_1 : \mu >\mu_1,
\] for some \(\mu_0 < \mu_0\) in the space of mean parameters. For any
constant \(g > 0\), Lorden 1973 proved the following upper bound on the
boundary-crossing probability:<sup>[1](#myfootnote1)</sup> \[
\sup_{\mu \leq \mu_0}\mathbb{P}_{\mu} \left(n \geq 1:   \log\mathrm{GLR}_n(>\mu_1,\leq\mu_0) \geq g\right)
\leq \begin{cases} e^{-g} &\mbox{if } \mathrm{KL}(\mu_1,\mu_0) \geq g \\
\left(1 + \frac{g}{\mathrm{KL}(\mu_1,\mu_0)} \right)e^{-g} & \mbox{otherwise }.
\end{cases}
\]

Alternatively, in our paper, we proved the following inequality holds:

\[
\sup_{\mu \leq \mu_0}\mathbb{P}_{\mu} \left(\exists n \geq 1:   \log\mathrm{GLR}_n(>\mu_1,\leq\mu_0) \geq g\right) 
\leq \begin{cases} e^{-g} &\mbox{if } \mathrm{KL}(\mu_1,\mu_0) \geq g \\
\inf_{\eta >1} \left\lceil \log_\eta \left(\frac{g}{\mathrm{KL}(\mu_1,\mu_0)}\right)\right\rceil e^{-g / \eta} & \mbox{otherwise }. 
\end{cases}
\] (In the paper, a nonparametric generalization of the above inequaltiy
was presented.)

Now, for any \(\alpha \in (0,1)\), let \(g_\alpha^L(\mu_1,\mu_0)\) and
\(g_\alpha(\mu_1,\mu_0)\) be smallest boundary values which make RHS of
two above inequalities equal to \(\alpha\), respectively. The following
Rcode reproduces the Fig. 3 in Section III-A in which we compared
\(g_\alpha^L(\mu_1,\mu_0)\) and \(g_\alpha(\mu_1,\mu_0)\) for normal
distributions with \(\sigma = 1\). (In this case,
\(\mathrm{KL}(\mu_1 ,\mu_0) = (\mu_1 - \mu_0)^2 / 2\).)

``` r
library(latex2exp)
library(SGLRT)


alpha <- 10^seq(-1,-10)
theta_vec <- c(1,0.5, 10^seq(-1,-10)) 
theta_vec <- 2^seq(0,-30)

d_vec <- theta_vec^2 / 2

for (i in seq_along(alpha[1:3])){
  f_lorden <- function(d) const_boundary_lorden(alpha[i], d)
  lorden <- sapply(d_vec, f_lorden)
  f_ours <- function(d) const_boundary(alpha[i], d)
  ours  <- sapply(d_vec, f_ours)
  if (i == 1){
    plot(1/theta_vec, unlist(lorden["g",]), type = "l", 
         log ="x",
         ylab = "Boundary Value",
         xlab = TeX("$|\\mu_1 - \\mu_0|^{-1}$ (log scale)"))
    points(1/theta_vec, unlist(ours["g",]), type = "l", col = 2)
  } else {
    points(1/theta_vec, unlist(lorden["g",]), type = "l", col = 1,
           lty = i)
    points(1/theta_vec, unlist(ours["g",]), type = "l", col = 2,
           lty = i)
  }
}
legend("topleft", TeX(c(paste0("Lorden's ($\\alpha = ", alpha[1:3],"$)"), paste0("Ours ($\\alpha = ", alpha[1:3],"$)"))), col = c(rep(1,3), rep(2,3)), lty = rep(1:3,2))
```

<img src="man/figures/README-Fig.3-1.png" width="100%" />

### 2\. Ratio of CI’s width to CLT *(Fig. 5 in Section IV-C)*

For each \(n\), let \(\mathrm{CI}_n(\alpha)\) be a random set based on
first \(n\) observations \(X_1, \dots, X_n\) from a distribution with
mean \(\mu\). If the sequence of random sets
\(\{\mathrm{CI}_n(\alpha)\}_{n \in \mathbb{N}}\) satisfies the following
inequality: \[
\mathbb{P}_\mu \left( \mu \in \mathrm{CI}_n(\alpha), \forall n \in \mathbb{N}\right) \geq 1-\alpha,~~\forall \mu \in (0,1),
\] it is called a *confidence sequence* with confidence level \(\alpha\)
(e.g, see [Howard et
al., 2018](https://arxiv.org/abs/1810.08240)<sup>[2](#myfootnote2)</sup>
and references therein for a comprehensive summary of related
literature.)

Based on the upper bound on the boundary crossing probability above, in
our paper, we presented two novel confidence sequences called *GLR-like*
and *discrete mixture* based confidence seqeunces. For sub-Gaussian
distributions with parameter \(\sigma\), The GLR-like confidence
seqeunce is given by \[
    \mathrm{CI}_n^{\mathrm{G}} := \begin{cases}  \left(\bar{X}_n - \sigma\sqrt{\frac{g_\alpha}{n_{\min}}}\left[  \sqrt{2}  + \frac{1}{\sqrt{2}}\left(\frac{n_{\min}}{n}-1\right)\right], \infty\right) &\mbox{if } n \in [n_0,n_{\min}) \\
 \left( \bar{X}_n - \sigma\sqrt{\frac{2g_\alpha}{n}}, \infty\right) &\mbox{if }  n \in [n_{\min}, n_{\max}] \\
 \left(\bar{X}_n - \sigma\sqrt{\frac{g_\alpha}{n_{\max}}}\left[ \sqrt{2}  + \frac{1}{\sqrt{2}}\left(\frac{n_{\max}}{n}-1\right)\right], \infty\right) &\mbox{if } n \in (n_{\max}, \infty)
  \end{cases},
\] for any given level \(\alpha \in (0, 1)\) and target time interval
\([n_{\min}, n_{\max}]\) on which the confidence seqeunce is
time-uniformly close to the Chrenoff bound. Unlike the GLR-like one, the
discrete mixture based confidence sequence
\(\{\mathrm{CI}_n^{\mathrm{DM}}\}\) does not have an explicit form and
see the paper for the detailed explanation how to compute
\(\mathrm{CI}_n^{\mathrm{DM}}\) for given sample mean \(\bar{X}_n\) and
confidence level \(\alpha\). By the construction, we have
\(\mathrm{CI}_n^{\mathrm{G}} \subset \mathrm{CI}_n^{\mathrm{DM}}\).

In Section IV-C of the paper, we compare these bounds with the stithcing
and normal mixture bounds in [Howard et
al., 2018](https://arxiv.org/abs/1810.08240)<sup>[2](#myfootnote2)</sup>
where each confidence intervals for stitching
\(\mathrm{CI}_{n}^{\mathrm{ST}}\) and normal mixture method
\(\mathrm{CI}_n^{\mathrm{NM}}\) is given by \[
\begin{aligned}
    \mathrm{CI}_{n}^{\mathrm{ST}} &:= \left(\bar{X}_n-  \frac{1.7}{\sqrt{n}}\sqrt{\log\log(2n) + 0.72 \log\left(\frac{5.2}{\alpha}\right)}, \infty \right) \\
    \mathrm{CI}_{n}^{\mathrm{NM}} &:= \left(\bar{X}_n-  \sqrt{2 \left(1 + \frac{\rho}{n}\right) \log\left(\frac{1}{2\alpha} \sqrt{\frac{n + \rho}{\rho} + 1}\right)}, \infty \right),
\end{aligned}
\] where we set \(\rho = 1260\) by following the setting in Figure 9 of
[Howard et
al., 2018](https://arxiv.org/abs/1810.08240)<sup>[2](#myfootnote2)</sup>.

The following Rcode reproduces the Fig. 5 in Section IV-C in which we
compared ratios of widths of confidence intervals above to the
poinstwise and asymptotically valid normal confidence intervals based on
the central limit theorem.

``` r
library(SGLRT)

# Hoeffding
hoeff <- function(v, alpha = 0.025){
  sqrt(2  * log(1/alpha) / v)
}

# Normal mixture function
normal_mix <- function(v, alpha = 0.025, rho = 1260){
  sqrt(2 * (1 + rho / v) * log(1/(2*alpha) * sqrt((v + rho)/rho) + 1)) 
}

#Stitching
stitch <- function(v, alpha = 0.025){
  1.7 * sqrt((log(log(2) + log(v)) + 0.72 * log(5.2/alpha)) / v) 
}

# CLT
CLT <- function(v, alpha = 0.025){
  qnorm(1-alpha) / sqrt(v)
}

# Functions to construct our bounds
# First one has target interval [1, nmax]
# Second one has target interval [nmax / 20, namx * 4]

alpha <- 0.025
nmax <- 1e+5

nmax1 <- nmax
nmin1 = 1

nmax2 <- nmax1 * 4
nmin2 = nmax1 / 20


ours1 <- SGLR_CI_additive(alpha, nmax1, nmin1)
ours2 <- SGLR_CI_additive(alpha, nmax2, nmin2)


# Compute the width of CIs on exponetially-spaced grid.
M <- log(nmax, base = 10)
v <- 10^(seq(0,M + 0.2, length.out = 100))


# Compute existing bounds
hoeff_vec <- sapply(v, hoeff)
CLT_vec <- sapply(v, CLT)
normal_mix_vec <- sapply(v, normal_mix)
stit_vec <- sapply(v, stitch)

existing_list <- list(stitch = stit_vec,
                      normal_mix = normal_mix_vec)

# Compute our bounds
GLR_like_1_vec <- sapply(v, ours1$GLR_like_fn)
GLR_like_2_vec <- sapply(v, ours2$GLR_like_fn)
our_dis_mix_1_vec <- sapply(v, ours1$dis_mix_fn)
our_dis_mix_2_vec <- sapply(v, ours2$dis_mix_fn)


ours_list_1 <- list(GLR_like_1 = GLR_like_1_vec,
                    our_dis_mix_1 = our_dis_mix_1_vec)

ours_list_2 <- list(GLR_like_2 = GLR_like_2_vec,
                    our_dis_mix_2 = our_dis_mix_2_vec)

# Plot ratio of bounds
 title <- "Ratio of CI's width to CLT"
 plot(v, hoeff_vec / CLT_vec, type = "l",
         main = title,
         ylab = "Ratio",
         xlab = "n",
         ylim = c(1, 4),
         xlim = c(1, nmax))
  col = 1
  legend_col <- c(1)
  for (i in seq_along(existing_list)){
    col <- col + 1
    legend_col <- c(legend_col, col)
    lines(v, existing_list[[i]] / CLT_vec, col = col)
  }
  legend_lty <- rep(1, length(existing_list) + 1)
   for (i in seq_along(ours_list_1)){
      col <- col + 1
      legend_col <- c(legend_col, col)
      lines(v, ours_list_1[[i]] / CLT_vec,
            lty = 2, lwd = 2, col = col)
    }
    legend_lty <- c(legend_lty, rep(2, length(ours_list_1)))
    for (i in seq_along(ours_list_2)){
      col <- col + 1
      legend_col <- c(legend_col, col)
      lines(v, ours_list_2[[i]] / CLT_vec,
            lty = 4, lwd = 2, col = col)
    }
    legend_lty <- c(legend_lty, rep(4, length(ours_list_2)))
    abline(v = c(nmin1, nmin2, nmax1, nmax2), lty = 3)
    bounds_name <- c("Hoeffiding", "Stitching", "Normal Mixture",
                "GLR-like 1", "Discrete Mixture 1",
                "GLR-like 2", "Discrete Mixture 2")
     legend("topright", bounds_name,
         lty = legend_lty,
         col = legend_col,
         bg= "white")
```

<img src="man/figures/README-Fig.5-1.png" width="100%" />

## Reference

<a name="myfootnote1">1</a>: G. Lorden, “Open-ended tests for
koopman-darmois families,”The Annals of Statistics, vol. 1, no. 4,
pp. 633–643, 1973.

<a name="myfootnote2">2</a>: S. R. Howard, A. Ramdas, J. McAuliffe, and
J. Sekhon, “Uniform, nonparametric, non-asymptotic confidence
sequences,”arXiv preprint arXiv:1810.08240, 2018.
