#' Fit flexible Cox models across multiple biomarkers
#'
#' @param data Input data frame.
#' @param biomarkers Character vector of biomarker base names.
#' @param event_col Event column name.
#' @param biomarker_type Biomarker mode.
#' @param time_col Time column for standard Cox.
#' @param start_col Start column for start-stop Cox.
#' @param stop_col Stop column for start-stop Cox.
#' @param covariates Optional covariate names.
#' @param interaction_var Optional interaction variable name.
#' @param time_varying_coefficients Whether to add a time-varying coefficient
#'   for each baseline biomarker. Supported only when `biomarker_type` is
#'   `"baseline"`.
#' @param center_time Center time if needed
#' @param baseline_suffix Baseline suffix.
#' @param tv_suffix Time-varying suffix.
#' @param change_suffix Change suffix.
#' @param ties Tie handling method.
#'
#' @return A list containing a stacked result data frame and term information.
#' @export

run_multiple_cox_flexible <- function(data,
                                      biomarkers,
                                      event_col,
                                      biomarker_type = c("baseline", "time_varying" ,"baseline_change"),
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
  results <- dplyr::bind_rows(
    lapply(biomarkers, function(biomarker) {
      tryCatch(
        run_single_cox_flexible(
          data = data,
          biomarker = biomarker,
          event_col = event_col,
          biomarker_type = biomarker_type,
          time_col = time_col,
          start_col = start_col,
          stop_col = stop_col,
          covariates = covariates,
          interaction_var = interaction_var,
	  time_varying_coefficients=time_varying_coefficients,
	  center_time=center_time,
          baseline_suffix = baseline_suffix,
          tv_suffix = tv_suffix,
          change_suffix = change_suffix,
          ties = ties
        ),
        error = function(e) {
          data.frame(
            biomarker = biomarker,
            term = NA_character_,
            coef = NA_real_,
            HR = NA_real_,
            lower95 = NA_real_,
            upper95 = NA_real_,
            pvalue = NA_real_,
            formula = NA_character_,
            effect_type = NA_character_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        }
      )
    })
  )

  rownames(results) <- NULL


  if (biomarker_type == "baseline") {
    main_terms <- c(baseline_suffix)
  }

  if (biomarker_type == "time_varying") {
    main_terms <- c(tv_suffix)
  }

  if (biomarker_type == "baseline_change") {
    main_terms <- c(baseline_suffix, change_suffix)
  }
 
  other_terms=c()
  if (!is.null(interaction_var)) {
    other_terms <- c(other_terms,paste0(main_terms, ":", interaction_var))
  }

  if (time_varying_coefficients){
    other_terms <- c(other_terms,main_terms)
  }
 term_info= list(
    main_terms = unique(main_terms),
    other_terms = unique(other_terms)
  )

  return(list(results, term_info))
}
