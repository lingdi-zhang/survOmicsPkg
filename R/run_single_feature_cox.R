#' Fit one flexible Cox model for one biomarker
#'
#' @param data Input data frame.
#' @param biomarker Biomarker base name.
#' @param event_col Event column name.
#' @param biomarker_type Biomarker mode.
#' @param time_col Time column for standard Cox.
#' @param start_col Start column for start-stop Cox.
#' @param stop_col Stop column for start-stop Cox.
#' @param covariates Optional covariate names.
#' @param interaction_var Optional interaction variable name.
#' @param time_varying_coefficients Whether to add a time-varying coefficient
#'   for a baseline biomarker.
#' @param center_time A finite numeric scalar when time-varying coefficients are enabled.
#' @param baseline_suffix Baseline suffix.
#' @param tv_suffix Time-varying suffix.
#' @param change_suffix Change suffix.
#' @param ties Tie handling method.
#'
#' @details The returned fdr_group column identifies baseline, change, interaction, and time coefficient families. Nonstandard column names are supported.
#' @details Longitudinal modes require both interval columns. Results include
#'   `fit_status` (ok, warning, unreliable, partial, or not_estimable) and
#'   captured `warning` messages.
#'   Convergence and infinite-coefficient warnings mark a fit as unreliable.
#'   Survival times must be finite, non-missing numeric values; events must be
#'   non-missing numeric 0/1 values. Interval models require stop > start.
#'   Interaction models include the modifier's main effect as a covariate.
#'   All biomarker, covariate, and modifier columns must exist in `data`.
#'   Baseline and change model columns must have distinct names.
#'   Each row has an `effect_status` of estimable or not_estimable, based on
#'   finite coefficients and p-values. Fits with some unestimable biomarker
#'   effects are partial; fits with none estimable are not_estimable.
#'   Convergence failures retain unreliable status. Covariate estimability
#'   does not determine the overall biomarker fit status.
#' @return A data frame of model output.
#' @export

run_single_cox_flexible <- function(data,
                                    biomarker,
                                    event_col,
                                    biomarker_type = c("baseline", "time_varying", "baseline_change"),
                                    time_col = NULL,
                                    start_col = NULL,
                                    stop_col = NULL,
                                    covariates = NULL,
                                    interaction_var = NULL,
				    time_varying_coefficients=FALSE,
				    center_time=0,
                                    baseline_suffix = "_bl",
                                    tv_suffix = "_tv",
                                    change_suffix = "_delta",
                                    ties = "efron") {
  biomarker_type <- match.arg(biomarker_type)

  if (!is.logical(time_varying_coefficients) ||
      length(time_varying_coefficients) != 1L ||
      is.na(time_varying_coefficients)) {
    stop("`time_varying_coefficients` must be TRUE or FALSE.", call. = FALSE)
  }
  if (time_varying_coefficients && biomarker_type != "baseline") {
    stop(
      "`time_varying_coefficients = TRUE` is supported only for baseline biomarkers.",
      call. = FALSE
    )
  }

  if (time_varying_coefficients &&
      (!is.numeric(center_time) || length(center_time) != 1L ||
       is.na(center_time) || !is.finite(center_time))) {
    stop("`center_time` must be one finite, non-missing numeric value.", call. = FALSE)
  }

  if (biomarker_type != "baseline" &&
      (is.null(start_col) || is.null(stop_col))) {
    stop("Longitudinal models require both `start_col` and `stop_col`.", call. = FALSE)
  }

  validate_survival_input(data, event_col, time_col, start_col, stop_col)

  biomarker_columns <- switch(biomarker_type,
    baseline = paste0(biomarker, baseline_suffix),
    time_varying = paste0(biomarker, tv_suffix),
    baseline_change = c(paste0(biomarker, baseline_suffix),
                        paste0(biomarker, change_suffix)))
  duplicate_columns <- unique(biomarker_columns[duplicated(biomarker_columns)])
  if (length(duplicate_columns)) {
    stop(paste0("Biomarker model columns must be distinct: ",
                paste(duplicate_columns, collapse = ", "), "."), call. = FALSE)
  }
  required_model_columns <- unique(c(biomarker_columns, covariates, interaction_var))
  missing_columns <- setdiff(required_model_columns, names(data))
  if (length(missing_columns)) {
    stop(paste0("Missing model columns in `data`: ",
                paste(missing_columns, collapse = ", "), "."), call. = FALSE)
  }

  term_info <- build_biomarker_terms(
    biomarker = biomarker,
    biomarker_type = biomarker_type,
    baseline_suffix = baseline_suffix,
    tv_suffix = tv_suffix,
    change_suffix = change_suffix,
    interaction_var = interaction_var
  )

  rhs_terms <- unique(c(term_info$rhs_terms, lapply(covariates, as.name)))

  fml <- build_surv_formula(
    event_col = event_col,
    time_col = time_col,
    start_col = start_col,
    stop_col = stop_col,
    rhs_terms = rhs_terms
  )

  if (time_varying_coefficients) {
    time_varying_term <- paste0(biomarker, baseline_suffix)
    fml[[3L]] <- call("+", fml[[3L]], call("tt", as.name(time_varying_term)))
  }
  fit_warnings <- character()
  fit <- withCallingHandlers(
    if (time_varying_coefficients) {
      survival::coxph(fml, data = data, ties = ties,
        tt = function(x, t, ...) x * (t - center_time))
    } else {
      survival::coxph(fml, data = data, ties = ties)
    },
    warning = function(w) {
      fit_warnings <<- c(fit_warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  unreliable <- any(grepl("converg|iterations|infinite", fit_warnings,
                          ignore.case = TRUE))
  fit_status <- if (unreliable) "unreliable" else
    if (length(fit_warnings)) "warning" else "ok"
  sm <- summary(fit)
  coefficient_model_terms <- rep(NA_character_, length(stats::coef(fit)))
  for (model_term in names(fit$assign)) {
    coefficient_model_terms[fit$assign[[model_term]]] <- model_term
  }
  names(coefficient_model_terms) <- names(stats::coef(fit))

  out <- data.frame(
    biomarker = biomarker,
    term = rownames(sm$coefficients),
    coef = sm$coefficients[, "coef"],
    HR = sm$conf.int[, "exp(coef)"],
    lower95 = sm$conf.int[, "lower .95"],
    upper95 = sm$conf.int[, "upper .95"],
    pvalue = sm$coefficients[, "Pr(>|z|)"],
    formula = paste(deparse(fml), collapse = ""),
    fit_status = fit_status,
    warning = if (length(fit_warnings)) paste(unique(fit_warnings), collapse = "; ") else NA_character_,
    row.names = NULL,
    stringsAsFactors = FALSE
  )

  model_terms <- unname(coefficient_model_terms[out$term])
  is_time_varying_term <- startsWith(model_terms, "tt(") &
    vapply(
      model_terms,
      function(term) {
        if (is.na(term) || !startsWith(term, "tt(")) {
          return(FALSE)
        }
        sub("^tt\\((.*)\\)$", "\\1", term) %in% term_info$main_terms
      },
      logical(1)
    )
  is_main_term <- model_terms %in% term_info$main_terms & !is_time_varying_term
  is_other_term <- model_terms %in% term_info$other_terms | is_time_varying_term
  out$effect_type <- ifelse(
    is_main_term, "main_terms",
    ifelse(is_other_term, "other_terms", "covariate")
  )
  main_groups <- switch(biomarker_type,
    baseline = "baseline", time_varying = "time_varying",
    baseline_change = c("baseline", "change"))
  group_by_term <- stats::setNames(main_groups, term_info$main_terms)
  if (length(term_info$other_terms)) {
    group_by_term <- c(group_by_term, stats::setNames(
      paste0("interaction_", main_groups), term_info$other_terms))
  }
  out$fdr_group <- unname(group_by_term[model_terms])
  out$fdr_group[is_time_varying_term] <- "time_coefficient"
  estimable <- is.finite(out$coef) & is.finite(out$pvalue)
  out$effect_status <- ifelse(estimable, "estimable", "not_estimable")
  requested_effects <- is_main_term | is_other_term
  # Convergence failures take precedence over missing individual effects.
  if (!unreliable && any(requested_effects & !estimable)) {
    out$fit_status <- if (any(requested_effects & estimable)) "partial" else "not_estimable"
  }
  out
}



# Validate direct fitting calls as well as preprocessed data.
validate_survival_input <- function(data, event_col, time_col, start_col, stop_col) {
  if (is.null(start_col) != is.null(stop_col)) {
    stop("Supply both `start_col` and `stop_col` together.", call. = FALSE)
  }
  interval <- !is.null(start_col)
  time_fields <- if (interval) c(start_col, stop_col) else time_col
  expected_time_fields <- if (interval) 2L else 1L
  fields <- c(time_fields, event_col)
  if (!is.data.frame(data) || length(event_col) != 1L ||
      length(time_fields) != expected_time_fields ||
      anyNA(fields) || any(!nzchar(fields)) || !all(fields %in% names(data))) {
    stop("Survival columns must name existing columns in `data`.", call. = FALSE)
  }
  events <- data[[event_col]]
  if (!is.numeric(events) || anyNA(events) || any(!is.finite(events)) ||
      any(!events %in% c(0, 1))) {
    stop("Event statuses must be non-missing numeric values coded as 0 or 1.", call. = FALSE)
  }
  for (field in time_fields) {
    times <- data[[field]]
    if (!is.numeric(times) || anyNA(times) || any(!is.finite(times))) {
      stop("Survival times must be non-missing, finite numeric values.", call. = FALSE)
    }
  }
  if (interval && any(data[[stop_col]] <= data[[start_col]])) {
    stop("Every survival interval must have stop > start.", call. = FALSE)
  }
}
