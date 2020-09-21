cat("evaluate binomial CIs")

library(tidyverse)

eval_dir <- "eval"

# Functions ===================================================================

prop_ci_cis <- function(success, total) {
  res <- PropCIs::exactci(success, total, 0.95)
  tibble(
    low = res$conf.int[[1]],
    point = success / total,
    high = res$conf.int[[2]],
    ci = "propci"
  )
}
beta_cis <- function(success, total) {
  tibble(
    low = qbeta(0.025, success + 1, total - success + 1),
    point = success / total,
    high = qbeta(0.975, success + 1, total - success + 1),
    ci = "beta"
  )
}
normal_cis <- function(success, total) {
  est <- success / total
  sd <- sqrt(est * (1 - est) / total)
  tibble(
    point = est,
    low = qnorm(0.025, point, sd),
    high = qnorm(0.975, point, sd),
    ci = "normal"
  )
}
cis <- function(success, total) {
  bind_rows(
    prop_ci_cis(success, total),
    beta_cis(success, total),
    normal_cis(success, total)
  )
}

# Script ======================================================================

p <- 0.5
t <- 20
N <- 1e4
tibble(
  sample = 1:N,
  success = rbinom(N, t, p),
  total = t,
) %>%
  group_by(sample, success, total) %>%
  summarise(cis(success, total), .groups = "drop") %>%
  group_by(ci) %>%
  summarise(sum(low < p & high > p) / n())
