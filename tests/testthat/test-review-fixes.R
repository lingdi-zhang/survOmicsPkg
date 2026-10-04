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
  expect_false("unused_feature" %in% names(result))
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

test_that("baseline event times use original statuses with default column names", {
  metadata <- data.frame(subject = c(1, 1, 1, 2, 2),
                         Visit = c(0, 2, 4, 0, 6), event = c(0, 0, 1, 0, 0))
  result <- construct_metadata(metadata, event_col = "event")
  expect_equal(result$time, c(4, 6))
  expect_equal(result$event, c(1L, 0L))
})

test_that("explicit FDR groups separate interactions and time coefficients", {
  set.seed(12)
  n <- 160
  data <- data.frame(time = rexp(n), event = rbinom(n, 1, .7),
                     marker_bl = rnorm(n), second_bl = rnorm(n),
                     group = factor(sample(c("a", "b"), n, TRUE)))
  fits <- run_multiple_cox_flexible(data, c("marker", "second", "missing"),
    event_col = "event", time_col = "time", interaction_var = "group",
    time_varying_coefficients = TRUE)
  result <- calculate_FDR(fits[[1]], fits[[2]])
  expect_equal(nrow(result), nrow(fits[[1]]))
  expect_equal(result$term, fits[[1]]$term)
  expect_setequal(na.omit(result$fdr_group),
                  c("baseline", "interaction_baseline", "time_coefficient"))
  for (group in c("baseline", "interaction_baseline", "time_coefficient")) {
    rows <- which(result$fdr_group == group)
    expect_length(rows, 2)
    expect_equal(result$FDR[rows], p.adjust(result$pvalue[rows], "BH"))
  }
  expect_true(is.na(result$FDR[result$biomarker == "missing"]))
})

test_that("baseline and change FDR groups do not depend on overlapping names", {
  set.seed(13)
  data <- data.frame(start = 0, stop = rexp(100), event = rbinom(100, 1, .7),
                     marker_bl = rnorm(100), marker_bl_delta = rnorm(100))
  result <- run_single_cox_flexible(data, "marker", "event",
    biomarker_type = "baseline_change", start_col = "start", stop_col = "stop",
    change_suffix = "_bl_delta")
  expect_equal(result$fdr_group, c("baseline", "change"))
  expect_equal(calculate_FDR(result)$FDR, result$pvalue)
})

test_that("longitudinal metadata rejects first-visit events and invalid times", {
  m <- data.frame(subject = 1, Visit = 0:2, event = c(1, 0, 0))
  for (mode in c("time_varying", "baseline_change")) {
    expect_error(construct_metadata(m, event_col = "event", biomarker_type = mode),
                 "First-visit events.*1")
  }
  m$event <- c(0, 0, 1)
  for (times in list(c(0, 0, 2), c(0, NA, 2), c(0, Inf, 2))) {
    m$Visit <- times
    expect_error(construct_metadata(m, event_col = "event", biomarker_type = "time_varying"),
                 "duplicate keys|finite numeric")
  }
})

test_that("preprocessing rejects duplicate feature and covariate keys", {
  m <- data.frame(subject = c(1, 1), Visit = c(0, 2), event = c(0, 1))
  b <- data.frame(subject = 1, marker = 2)
  tv <- data.frame(subject = c(1, 1), Visit = c(0, 2), marker = c(2, 3))
  expect_error(preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = rbind(b, b)), "baseline_feature_table.*duplicate keys")
  expect_error(preprocess_data(m, "marker", event_col = "event",
    biomarker_type = "time_varying", time_varying_feature_table = rbind(tv, tv)),
    "time_varying_feature_table.*duplicate keys")
  expect_error(preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b, covariates_table = data.frame(subject = c(1, 1), age = c(20, 21))),
    "covariates_table.*duplicate keys")
})

test_that("nonstandard names work in survival, interaction, and time coefficient terms", {
  set.seed(20)
  n <- 180
  d <- data.frame(check.names = FALSE, "visit time" = rexp(n),
    "event status" = rbinom(n, 1, .7), "IL-6_bl" = rnorm(n),
    "patient group" = factor(sample(c("a", "b"), n, TRUE)), "patient age" = rnorm(n))
  result <- run_single_cox_flexible(d, "IL-6", "event status", time_col = "visit time",
    interaction_var = "patient group", covariates = "patient age", time_varying_coefficients = TRUE)
  expect_setequal(na.omit(result$fdr_group), c("baseline", "interaction_baseline", "time_coefficient"))
  expect_equal(sum(result$effect_type == "main_terms"), 1)
  expect_equal(sum(result$effect_type == "other_terms"), 2)
  expect_true(all(is.finite(result$pvalue)))
  # Compare estimates against a directly specified survival model.
  direct <- survival::coxph(survival::Surv(`visit time`, `event status`) ~
    `IL-6_bl` + `patient group` + `IL-6_bl`:`patient group` + `patient age` + tt(`IL-6_bl`),
    data = d, tt = function(x, t, ...) x * t)
  expect_equal(result$coef, unname(coef(direct)))
})

test_that("legacy FDR uses structural suffixes and rejects ambiguity", {
  r <- data.frame(term = c("marker_delta_bl", "other_bl", "marker_delta_delta", "other_delta"),
    effect_type = "main_terms", pvalue = c(.01, .03, .8, .9))
  terms <- list(main_terms = c("_bl", "_delta"), other_terms = character())
  result <- calculate_FDR(r, terms)
  expect_equal(result$FDR, c(.02, .03, .9, .9))
  expect_equal(nrow(result), nrow(r))
  terms$main_terms <- c("_bl", "bl")
  expect_error(calculate_FDR(r, terms), "Ambiguous FDR group")
  r <- data.frame(term = c("tt(marker_delta_bl)", "marker_delta_bl:groupb"),
    effect_type = "other_terms", pvalue = c(.1, .2))
  result <- calculate_FDR(r, list(main_terms = c("_bl", "_delta"),
    other_terms = c("_bl:group", "_delta:group", "_bl", "_delta")))
  expect_equal(result$fdr_group, c("time:_bl", "other:_bl:group"))
})


test_that("incomplete interval inputs never fall back to standard survival", {
  d <- data.frame(start = 0, stop = 1:4, event = c(0, 1, 0, 1), marker_tv = 1:4)
  for (mode in c("time_varying", "baseline_change")) {
    for (cols in list(list(), list(start_col = "start"), list(stop_col = "stop"))) {
      expect_error(do.call(run_single_cox_flexible,
        c(list(data = d, biomarker = "marker", event_col = "event",
               biomarker_type = mode, time_col = "stop"), cols)),
        "require both")
    }
  }
  expect_error(build_surv_formula("event", time_col = "stop",
    start_col = "start", rhs_terms = "marker_tv"), "both")
  expect_error(build_surv_formula("event", time_col = "stop",
    stop_col = "stop", rhs_terms = "marker_tv"), "both")
})

test_that("unused feature columns cannot collide with model fields", {
  m <- data.frame(subject = c(1, 1, 2, 2), Visit = c(0, 1, 0, 1), event = c(0, 1, 0, 0))
  b <- data.frame(subject = 1:2, marker = 2:3, time = 0, event = 0, age = 99)
  tv <- data.frame(subject = m$subject, Visit = m$Visit, marker = 1:4,
                   event = 0, stop = 99, age = 99)
  cv <- data.frame(subject = 1:2, age = c(20, 30))
  for (mode in c("baseline", "time_varying", "baseline_change")) {
    out <- preprocess_data(m, "marker", event_col = "event",
      baseline_feature_table = b, time_varying_feature_table = tv,
      covariates_table = cv, biomarker_type = mode)
    expect_equal(out$age, c(20, 30))
    expect_equal(out$event, c(1L, 0L))
    expect_false(any(grepl("\\.[xy]$", names(out))))
  }
  expect_error(preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b,
    covariates_table = data.frame(subject = 1:2, event = 0)),
    "Overlapping non-key columns: event")
})

test_that("nonconvergence is retained and excluded from FDR", {
  d <- data.frame(time = c(6, 6), event = c(1, 0), marker_bl = c(2, 3))
  bad <- run_single_cox_flexible(d, "marker", "event", time_col = "time")
  expect_equal(bad$fit_status, "unreliable")
  expect_match(bad$warning, "converg|iterations|infinite")
  expect_true(is.finite(bad$pvalue))
  good <- bad
  good$fit_status <- "ok"
  good$warning <- NA_character_
  good$pvalue <- .02
  result <- calculate_FDR(rbind(bad, good))
  expect_true(is.na(result$FDR[1]))
  expect_equal(result$FDR[2], .02)
  batch <- run_multiple_cox_flexible(d, c("marker", "missing"),
    "event", time_col = "time")[[1]]
  expect_equal(batch$fit_status, c("unreliable", "error"))
  expect_true(all(is.na(calculate_FDR(batch)$FDR)))
})


test_that("direct fits reject invalid survival fields and batch fits record errors", {
  set.seed(30)
  d <- data.frame(start = 0, stop = rexp(100) + 1,
    event = rbinom(100, 1, .7), marker_tv = rnorm(100))
  for (bad_event in list(NA_real_, Inf, 2, "1", TRUE)) {
    bad <- d
    bad$event <- rep(bad_event, nrow(d))
    expect_error(run_single_cox_flexible(bad, "marker", "event",
      time_col = "stop"), "Event statuses")
  }
  for (field in c("start", "stop")) {
    for (bad_time in list(NA_real_, Inf, "1")) {
      bad <- d
      bad[[field]] <- rep(bad_time, nrow(d))
      expect_error(run_single_cox_flexible(bad, "marker", "event",
        biomarker_type = "time_varying", start_col = "start", stop_col = "stop"),
        "Survival times")
    }
  }
  expect_error(run_single_cox_flexible(transform(d, stop = NA_real_),
    "marker", "event", time_col = "stop"), "Survival times")
  expect_error(run_single_cox_flexible(d, "marker", "missing", time_col = "stop"),
    "existing columns")
  for (offset in c(0, 1)) {
    bad <- d
    bad$start[1] <- bad$stop[1] + offset
    expect_error(run_single_cox_flexible(bad, "marker", "event",
      biomarker_type = "time_varying", start_col = "start", stop_col = "stop"),
      "stop > start")
    batch <- run_multiple_cox_flexible(bad, c("marker", "second"), "event",
      biomarker_type = "time_varying", start_col = "start", stop_col = "stop")[[1]]
    expect_true(all(batch$fit_status == "error"))
    expect_true(all(grepl("stop > start", batch$error)))
    expect_true(all(is.na(calculate_FDR(batch)$FDR)))
  }
})

test_that("interaction main effects agree with direct numeric and factor models", {
  set.seed(31)
  d <- data.frame(time = rexp(200), event = rbinom(200, 1, .7),
    marker_bl = rnorm(200), marker_delta = rnorm(200), start = 0)
  for (modifier in list(rnorm(200), factor(rep(c("a", "b", "c", "d"), 50)))) {
    d$group <- modifier
    for (mode in c("baseline", "baseline_change")) {
      args <- if (mode == "baseline") list(time_col = "time") else
        list(start_col = "start", stop_col = "time")
      result <- do.call(run_single_cox_flexible, c(list(data = d,
        biomarker = "marker", event_col = "event", biomarker_type = mode,
        interaction_var = "group"), args))
      duplicated <- do.call(run_single_cox_flexible, c(list(data = d,
        biomarker = "marker", event_col = "event", biomarker_type = mode,
        interaction_var = "group", covariates = "group"), args))
      direct <- if (mode == "baseline") {
        survival::coxph(survival::Surv(time, event) ~ marker_bl * group, data = d)
      } else {
        survival::coxph(survival::Surv(start, time, event) ~
          (marker_bl + marker_delta) * group, data = d)
      }
      expect_equal(result$coef, unname(coef(direct)[result$term]))
      expect_equal(duplicated, result)
      group_rows <- startsWith(result$term, "group")
      expect_true(any(group_rows))
      expect_true(all(result$effect_type[group_rows] == "covariate"))
      expect_true(all(is.na(result$fdr_group[group_rows])))
    }
  }
})

test_that("metadata output names preserve original statuses and subject IDs", {
  m <- data.frame(subject = c(1, 1, 1, 2, 2), Visit = c(0, 2, 4, 0, 6),
    status = c(0, 0, 1, 0, 0))
  result <- construct_metadata(m, event_col = "status",
    new_time_col = "status", new_event_col = "Visit")
  expect_equal(result$status, c(4, 6))
  expect_equal(result$Visit, c(1L, 0L))
  for (args in list(list(new_time_col = "subject"),
                   list(new_time_col = "same", new_event_col = "same"),
                   list(new_time_col = ""), list(new_event_col = NULL))) {
    expect_error(do.call(construct_metadata, c(list(metadata = m,
      event_col = "status"), args)), "Output column names")
  }
  tv <- construct_metadata(m, event_col = "status", biomarker_type = "time_varying",
    start_col = "Visit", stop_col = "next")
  expect_equal(tv$Visit, c(0, 2, 0))
  expect_equal(tv[["next"]], c(2, 4, 6))
  expect_equal(tv$status, c(0L, 1L, 0L))
  expect_error(construct_metadata(m, event_col = "status", biomarker_type = "time_varying",
    start_col = "status"), "distinct")
  # Internal names in the caller's table also remain safe.
  names(m) <- c(".subject", ".visit", ".status")
  internal <- construct_metadata(m, id_col = ".subject", time_col = ".visit",
    event_col = ".status", new_time_col = ".status", new_event_col = ".visit")
  expect_equal(internal$.status, c(4, 6))
  expect_equal(internal$.visit, c(1L, 0L))
})

test_that("model inputs must come from data even when globals match", {
  check_globals <- function() {
    external <- c("review_external_marker_bl", "review_external_marker_tv",
                  "review_external_marker_delta", "review_external_age", "review_external_group")
    existed <- vapply(external, exists, logical(1), envir = .GlobalEnv, inherits = FALSE)
    previous <- mget(external[existed], envir = .GlobalEnv)
    on.exit({
      for (name in external[existed]) assign(name, previous[[name]], envir = .GlobalEnv)
      rm(list = external[!existed], envir = .GlobalEnv)
    })
    set.seed(40)
    d <- data.frame(time = rexp(160) + 1, start = 0,
                    event = rbinom(160, 1, .7), marker_bl = rnorm(160))
    for (name in external) assign(name, rnorm(nrow(d)), envir = .GlobalEnv)
    for (mode in c("baseline", "time_varying", "baseline_change")) {
      args <- if (mode == "baseline") list(time_col = "time") else
        list(start_col = "start", stop_col = "time")
      expect_error(do.call(run_single_cox_flexible, c(list(data = d,
        biomarker = "review_external_marker", event_col = "event", biomarker_type = mode), args)),
        "Missing model columns.*review_external_marker")
      batch <- do.call(run_multiple_cox_flexible, c(list(data = d,
        biomarkers = "review_external_marker", event_col = "event", biomarker_type = mode), args))[[1]]
      expect_equal(batch$fit_status, "error")
      expect_match(batch$error, "Missing model columns")
      expect_true(is.na(calculate_FDR(batch)$FDR))
    }
    expect_error(run_single_cox_flexible(d, "marker", "event", time_col = "time",
      covariates = "review_external_age"), "Missing model columns.*review_external_age")
    expect_error(run_single_cox_flexible(d, "marker", "event", time_col = "time",
      interaction_var = "review_external_group"), "Missing model columns.*review_external_group")
    # A mixed batch still preserves successful fits.
    batch <- run_multiple_cox_flexible(d, c("marker", "review_external_marker"),
      "event", time_col = "time")[[1]]
    expect_equal(batch$fit_status, c("ok", "error"))
  }
  check_globals()
})

test_that("generated feature names cannot overwrite model or join fields", {
  m <- data.frame(subject = c(1, 1, 2, 2), Visit = c(0, 2, 0, 2),
                   event = c(0, 1, 0, 0))
  b <- data.frame(subject = 1:2, event = c(2, 3), marker = c(2, 3))
  tv <- data.frame(subject = 1:2, Visit = 0, event = c(2, 4), marker = c(2, 4))
  expect_error(preprocess_data(m, "event", event_col = "event",
    baseline_feature_table = b, time_varying_feature_table = tv,
    biomarker_type = "baseline_change", change_suffix = ""),
    "Generated feature column names conflict: event")
  for (suffix in c("_bl", "_tv")) {
    expect_error(preprocess_data(m, "marker", event_col = "event",
      baseline_feature_table = b, time_varying_feature_table = tv,
      biomarker_type = "baseline_change", change_suffix = suffix),
      "Generated feature column names conflict: marker")
  }
  for (reserved in c("subject", "start", "stop", "Visit")) {
    expect_error(preprocess_data(m, "marker", event_col = "event",
      baseline_feature_table = b, time_varying_feature_table = tv,
      biomarker_type = "baseline_change", change_suffix = "",
      # Make each reserved field equal the generated name in turn.
      id_col = if (reserved == "subject") "marker" else "subject",
      start_col = if (reserved == "start") "marker" else "start",
      stop_col = if (reserved == "stop") "marker" else "stop",
      time_col = if (reserved == "Visit") "marker" else "Visit"),
      "Generated feature column names conflict: marker")
  }
  expect_error(preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b, baseline_suffix = "", new_time_col = "marker"),
    "Generated feature column names conflict: marker")
  expect_error(preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b, covariates_table = data.frame(subject = 1:2, marker_bl = 1)),
    "Generated feature column names conflict: marker_bl")
  # An empty suffix remains valid when the generated name is unused.
  out <- preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b, time_varying_feature_table = tv,
    biomarker_type = "baseline_change", change_suffix = "")
  expect_equal(out$marker, c(0, 1))
  expect_equal(out$event, c(1L, 0L))
  expect_equal(out$marker_bl, c(2, 3))
  baseline <- preprocess_data(m, "marker", event_col = "event",
    baseline_feature_table = b, baseline_suffix = "")
  expect_equal(baseline$marker, c(2, 3))
})

test_that("direct baseline-change models reject duplicate column names", {
  set.seed(70)
  d <- data.frame(start = 0, stop = rexp(160) + 1,
                   event = rbinom(160, 1, .7), marker_bl = rnorm(160))
  expect_error(run_single_cox_flexible(d, "marker", "event",
    biomarker_type = "baseline_change", start_col = "start", stop_col = "stop",
    change_suffix = "_bl"), "Biomarker model columns must be distinct: marker_bl")
  batch <- run_multiple_cox_flexible(d, "marker", "event",
    biomarker_type = "baseline_change", start_col = "start", stop_col = "stop",
    change_suffix = "_bl")[[1]]
  expect_equal(batch$fit_status, "error")
  expect_equal(batch$effect_status, "not_estimable")
  expect_match(batch$error, "must be distinct")
  expect_true(is.na(calculate_FDR(batch)$FDR))
})

test_that("constant biomarkers are explicitly not estimable", {
  set.seed(71)
  d <- data.frame(time = rexp(160), event = rbinom(160, 1, .7),
                   marker_bl = 1, age = rnorm(160))
  result <- run_single_cox_flexible(d, "marker", "event", time_col = "time",
                                     covariates = "age")
  marker <- result$term == "marker_bl"
  expect_true(all(result$fit_status == "not_estimable"))
  expect_equal(result$effect_status[marker], "not_estimable")
  expect_equal(result$effect_status[!marker], "estimable")
  expect_true(is.na(result$pvalue[marker]))
  expect_true(all(is.na(calculate_FDR(result)$FDR)))
  # An unestimable covariate does not invalidate an estimable biomarker.
  d$marker_bl <- rnorm(nrow(d))
  d$age <- 1
  result <- run_single_cox_flexible(d, "marker", "event", time_col = "time",
                                     covariates = "age")
  expect_true(all(result$fit_status == "ok"))
  expect_equal(result$effect_status[result$term == "age"], "not_estimable")
  expect_true(is.finite(calculate_FDR(result)$FDR[result$term == "marker_bl"]))
})

test_that("partial models preserve estimable effects and FDR families", {
  set.seed(72)
  d <- data.frame(start = 0, stop = rexp(200) + 1, event = rbinom(200, 1, .7),
                   marker_bl = rnorm(200), second_bl = rnorm(200))
  d$marker_delta <- d$marker_bl
  d$second_delta <- rnorm(nrow(d))
  batch <- run_multiple_cox_flexible(d, c("marker", "second"), "event",
    biomarker_type = "baseline_change", start_col = "start", stop_col = "stop")[[1]]
  expect_true(all(batch$fit_status[batch$biomarker == "marker"] == "partial"))
  expect_true(all(batch$fit_status[batch$biomarker == "second"] == "ok"))
  missing <- batch$term == "marker_delta"
  expect_equal(batch$effect_status[missing], "not_estimable")
  result <- calculate_FDR(batch)
  expect_true(is.na(result$FDR[missing]))
  baseline <- result$fdr_group == "baseline"
  expect_equal(result$FDR[baseline], p.adjust(result$pvalue[baseline], "BH"))
  expect_equal(result$FDR[result$term == "second_delta"],
               result$pvalue[result$term == "second_delta"])
  direct <- survival::coxph(survival::Surv(start, stop, event) ~
                             marker_bl + marker_delta, data = d)
  expect_equal(result$coef[result$biomarker == "marker"], unname(coef(direct)))
  # Convergence warnings still exclude even finite effects.
  bad <- run_single_cox_flexible(data.frame(time = c(6, 6), event = c(1, 0),
    marker_bl = c(2, 3)), "marker", "event", time_col = "time")
  expect_equal(bad$fit_status, "unreliable")
  expect_true(all(is.na(calculate_FDR(bad)$FDR)))
})

test_that("duplicate batch biomarkers fail before fitting", {
  # Missing data would fail a fit; duplicates must be rejected before that.
  expect_error(run_multiple_cox_flexible(NULL, c("a", "b", "a"), "event"),
    "Duplicate biomarker names: a")
  expect_error(run_multiple_cox_flexible(NULL, c("a", "b", "a", "b"), "event"),
    "Duplicate biomarker names: a, b")
})

test_that("time coefficients require a finite numeric scalar center", {
  set.seed(80)
  d <- data.frame(time = rexp(180), event = rbinom(180, 1, .7),
                   marker_bl = rnorm(180), second_bl = rnorm(180))
  for (center in list(c(0, 10), NULL, NA_real_, Inf, -Inf, NaN, "0", TRUE)) {
    expect_error(run_single_cox_flexible(d, "marker", "event", time_col = "time",
      time_varying_coefficients = TRUE, center_time = center),
      "center_time.*one finite")
    batch <- run_multiple_cox_flexible(d, c("marker", "second"), "event",
      time_col = "time", time_varying_coefficients = TRUE, center_time = center)[[1]]
    expect_true(all(batch$fit_status == "error"))
    expect_true(all(grepl("center_time.*one finite", batch$error)))
    expect_true(all(is.na(calculate_FDR(batch)$FDR)))
  }
  for (center in c(-.5, 0, .25)) {
    result <- run_single_cox_flexible(d, "marker", "event", time_col = "time",
      time_varying_coefficients = TRUE, center_time = center)
    direct <- survival::coxph(survival::Surv(time, event) ~ marker_bl + tt(marker_bl),
      data = d, tt = function(x, t, ...) x * (t - center))
    expect_equal(result$coef, unname(coef(direct)))
    expect_equal(result$pvalue, unname(summary(direct)$coefficients[, "Pr(>|z|)"]))
    expect_true(all(result$fit_status == "ok"))
  }
  # An unused center does not affect a standard model.
  expect_equal(run_single_cox_flexible(d, "marker", "event", time_col = "time",
    center_time = NULL), run_single_cox_flexible(d, "marker", "event", time_col = "time"))
})

test_that("FDR rejects invalid p-values before correcting any family", {
  base <- data.frame(fdr_group = "baseline", pvalue = c(.01, .1), fit_status = "ok")
  for (value in c(-.1, 1.1, Inf, -Inf)) {
    bad <- base
    bad$pvalue[2] <- value
    expect_error(calculate_FDR(bad), "finite values in")
    # Ineligible rows still must not contain invalid observed p-values.
    bad$fit_status[2] <- "error"
    expect_error(calculate_FDR(bad), "finite values in")
  }
  for (values in list(c("0.01", "0.1"), factor(c("0.01", "0.1")), c(TRUE, FALSE))) {
    bad <- base
    bad$pvalue <- values
    expect_error(calculate_FDR(bad), "must be numeric")
  }
  expect_error(calculate_FDR(base, pvalue_col = "missing"), "existing column")
  expect_error(calculate_FDR(base, pvalue_col = c("pvalue", "missing")), "existing column")
  expect_error(calculate_FDR(base, pvalue_col = NA_character_), "existing column")
})

test_that("FDR accepts boundaries, missing values, and custom p-value columns", {
  r <- data.frame(fdr_group = "baseline", p = c(0, .01, 1, NA_real_, NaN))
  result <- calculate_FDR(r, pvalue_col = "p")
  expect_equal(result$FDR, p.adjust(r$p, method = "BH"))
  expect_equal(nrow(result), nrow(r))
  expect_true(all(is.na(result$FDR[4:5])))
  empty <- r[FALSE, ]
  expect_equal(nrow(calculate_FDR(empty, pvalue_col = "p")), 0L)
})
