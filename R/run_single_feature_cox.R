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
#' @param center_time Center time if needed
#' @param baseline_suffix Baseline suffix.
#' @param tv_suffix Time-varying suffix.
#' @param change_suffix Change suffix.
#' @param ties Tie handling method.
#'
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

  term_info <- build_biomarker_terms(
    biomarker = biomarker,
    biomarker_type = biomarker_type,
    baseline_suffix = baseline_suffix,
    tv_suffix = tv_suffix,
    change_suffix = change_suffix,
    interaction_var = interaction_var
  )

  rhs_terms <- unique(c(term_info$rhs_terms, covariates))

  fml <- build_surv_formula(
    event_col = event_col,
    time_col = time_col,
    start_col = start_col,
    stop_col = stop_col,
    rhs_terms = rhs_terms
  )

  if (time_varying_coefficients){
	  time_varying_term <- paste0(biomarker, baseline_suffix)
	  fml <- stats::as.formula(paste0(
	    paste(deparse(fml), collapse = ""),
	    " + tt(", time_varying_term, ")"
	  ))
	  fit <- survival::coxph(
          formula = fml,
          data = data,
	  tt = function(x, t, ...) x * (t - center_time),
          ties = ties
  )
  }
  else{
  fit <- survival::coxph(
    formula = fml,
    data = data,
    ties = ties
  )
  }
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
 out
}


