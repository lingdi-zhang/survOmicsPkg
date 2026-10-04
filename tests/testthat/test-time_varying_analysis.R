test_that("longitudinal workflow agrees with a directly specified Cox model", {
  set.seed(102)
  n <- 180
  followup <- rexp(n) + 1
  event <- rbinom(n, 1, .7)
  metadata <- data.frame(subject = rep(seq_len(n), each = 3),
    time = as.vector(rbind(0, followup / 2, followup)),
    event = as.vector(rbind(0, 0, event)))
  features <- data.frame(subject = metadata$subject, time = metadata$time,
                         biomarker = rnorm(3 * n))
  out <- preprocess_data(metadata, "biomarker", event_col = "event",
    time_col = "time", time_varying_feature_table = features,
    biomarker_type = "time_varying")
  expect_equal(nrow(out), 2 * n)
  expect_equal(sum(out$event), sum(event))
  res <- run_multiple_cox_flexible(out, "biomarker", "event",
    biomarker_type = "time_varying", start_col = "start", stop_col = "stop")
  direct <- survival::coxph(survival::Surv(start, stop, event) ~ biomarker_tv, data = out)
  expect_equal(res[[1]]$coef, unname(coef(direct)))
  expect_equal(res[[1]]$pvalue, unname(summary(direct)$coefficients[, "Pr(>|z|)"]))
  expect_equal(res[[1]]$fit_status, "ok")
  expect_true(is.na(res[[1]]$warning))
  expect_equal(calculate_FDR(res[[1]])$FDR, res[[1]]$pvalue)
})
