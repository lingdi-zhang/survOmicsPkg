test_that("preprocessing ignores missing values in unused feature columns", {
  metadata <- data.frame(
    subject = c(1, 1, 2, 2, 3, 3),
    time = c(0, 1, 0, 1, 0, 1),
    event = c(0, 0, 0, 1, 0, 0)
  )
  baseline <- data.frame(
    subject = c(1, 2, 3),
    marker = c(2, 3, NA_real_),
    unused_feature = c(NA_real_, 4, 5)
  )

  result <- preprocess_data(
    metadata,
    biomarkers = "marker",
    event_col = "event",
    time_col = "time",
    baseline_feature_table = baseline,
    biomarker_type = "baseline"
  )

  expect_equal(result$subject, c(1, 2))
  expect_true(is.na(result$unused_feature[1]))
})

test_that("factor-expanded interaction coefficients retain their effect type", {
  set.seed(4)
  n <- 160
  data <- data.frame(
    time = rexp(n),
    event = rbinom(n, 1, 0.7),
    marker_bl = rnorm(n),
    group = factor(sample(c("a", "b"), n, replace = TRUE))
  )

  result <- run_single_cox_flexible(
    data,
    biomarker = "marker",
    event_col = "event",
    time_col = "time",
    interaction_var = "group"
  )

  interaction <- grepl("^marker_bl:group", result$term)
  expect_true(any(interaction))
  expect_true(all(result$effect_type[interaction] == "other_terms"))

  fdr_result <- calculate_FDR(
    result,
    list(main_terms = "_bl", other_terms = "_bl:group")
  )
  expect_true(any(grepl("^marker_bl:group", fdr_result$term)))
})

test_that("time-varying coefficients are restricted to baseline biomarkers", {
  data <- data.frame(
    start = c(0, 1),
    stop = c(1, 2),
    event = c(0, 1),
    marker_tv = c(1, 2)
  )

  expect_error(
    run_single_cox_flexible(
      data,
      biomarker = "marker",
      event_col = "event",
      biomarker_type = "time_varying",
      start_col = "start",
      stop_col = "stop",
      time_varying_coefficients = TRUE
    ),
    "supported only for baseline biomarkers"
  )
})

test_that("baseline time-varying coefficient is identified as an additional effect", {
  set.seed(4)
  n <- 160
  data <- data.frame(
    time = rexp(n),
    event = rbinom(n, 1, 0.7),
    marker_bl = rnorm(n)
  )

  result <- run_single_cox_flexible(
    data,
    biomarker = "marker",
    event_col = "event",
    time_col = "time",
    time_varying_coefficients = TRUE
  )

  expect_equal(
    result$effect_type[result$term == "tt(marker_bl)"],
    "other_terms"
  )
})

test_that("multiple biomarker results combine successful and failed fits", {
  set.seed(5)
  n <- 160
  data <- data.frame(
    time = rexp(n),
    event = rbinom(n, 1, 0.7),
    marker_bl = rnorm(n)
  )

  result <- run_multiple_cox_flexible(
    data,
    biomarkers = c("marker", "missing"),
    event_col = "event",
    time_col = "time",
    biomarker_type = "baseline"
  )[[1]]

  expect_true(any(result$biomarker == "marker" & !is.na(result$pvalue)))
  expect_true(any(result$biomarker == "missing" & !is.na(result$error)))
})

test_that("metadata construction rejects missing or invalid event statuses", {
  metadata <- data.frame(
    subject = c(1, 1, 2, 2),
    time = c(0, 1, 0, 1),
    event = c(0, 1, 0, 0)
  )

  for (biomarker_type in c("baseline", "time_varying")) {
    expect_error(
      construct_metadata(
        transform(metadata, event = c(0, NA, 0, 0)),
        event_col = "event",
        time_col = "time",
        biomarker_type = biomarker_type
      ),
      "non-missing numeric values coded as 0 or 1"
    )
    expect_error(
      construct_metadata(
        transform(metadata, event = c(0, 2, 0, 0)),
        event_col = "event",
        time_col = "time",
        biomarker_type = biomarker_type
      ),
      "non-missing numeric values coded as 0 or 1"
    )
  }
})

test_that("valid event statuses produce time-varying intervals", {
  metadata <- data.frame(
    subject = c(1, 1, 1),
    time = c(0, 1, 2),
    event = c(0, 0, 1)
  )

  result <- construct_metadata(
    metadata,
    event_col = "event",
    time_col = "time",
    biomarker_type = "time_varying"
  )

  expect_equal(result$event, c(0L, 1L))
  expect_equal(result$start, c(0, 1))
  expect_equal(result$stop, c(1, 2))
})
