test_that("baseline workflow agrees with a directly specified Cox model", {
  set.seed(101)
  n <- 180
  followup <- rexp(n)
  event <- rbinom(n, 1, .7)
  metadata <- data.frame(subject = rep(seq_len(n), each = 2),
    time = as.vector(rbind(0, followup)), event = as.vector(rbind(0, event)))
  baseline <- data.frame(subject = seq_len(n), biomarker = rnorm(n))
  out <- preprocess_data(metadata, "biomarker", event_col = "event",
    time_col = "time", baseline_feature_table = baseline)
  expect_equal(nrow(out), n)
  res <- run_multiple_cox_flexible(out, "biomarker", "event", time_col = "time")
  direct <- survival::coxph(survival::Surv(time, event) ~ biomarker_bl, data = out)
  expect_equal(res[[1]]$coef, unname(coef(direct)))
  expect_equal(res[[1]]$pvalue, unname(summary(direct)$coefficients[, "Pr(>|z|)"]))
  expect_equal(res[[1]]$fit_status, "ok")
  expect_true(is.na(res[[1]]$warning))
  expect_equal(calculate_FDR(res[[1]])$FDR, res[[1]]$pvalue)
})
