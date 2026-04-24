#' Build biomarker terms for  Cox models
#'
#' @param biomarker Biomarker base name.
#' @param biomarker_type One of "baseline", "time_varying",or "baseline_change".
#' @param baseline_suffix Suffix for baseline biomarker columns.
#' @param tv_suffix Suffix for time-varying biomarker columns.
#' @param change_suffix Suffix for change biomarker columns.
#' @param interaction_var Optional interaction variable name.
#'
#' @return A list with main_terms, other_terms, and rhs_terms.
#' @export

build_biomarker_terms <- function(biomarker,
                                  biomarker_type = c("baseline", "time_varying", "baseline_change"),
                                  baseline_suffix = "_bl",
                                  tv_suffix = "_tv",
                                  change_suffix = "_delta",
                                  interaction_var = NULL
				 # time_varying_coefficients=FALSE,
				 # center_time=0
				  ) {

  biomarker_type <- match.arg(biomarker_type)
  bl <- paste0(biomarker, baseline_suffix)
  tv <- paste0(biomarker, tv_suffix)
  ch <- paste0(biomarker, change_suffix)

  main_terms <- character(0)
  other_terms <- character(0)

  if (biomarker_type == "baseline") {
    main_terms <- bl
  }

  if (biomarker_type == "time_varying") {
    main_terms <- tv
  }

  if (biomarker_type == "baseline_change") {
    main_terms <- c(bl, ch)
  }
  if (!is.null(interaction_var)) {
    other_terms <- c(other_terms, paste0(main_terms, ":", interaction_var))
  }

#  if (time_varying_coefficients){
#    other_terms<-c(other_terms, paste0("tt(",main_terms,") ,tt = function(x, t, ...) x * (t -", center_time,")"))
#  }

 out= list(
    main_terms = unique(main_terms),
    other_terms = unique(other_terms),
    rhs_terms = unique(c(main_terms, other_terms))
  )

}



