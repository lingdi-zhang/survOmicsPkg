#' Build survival formula
#'
#' @param event_col Event column name.
#' @param time_col Time column name for standard Cox model.
#' @param start_col Start column name for start-stop Cox model.
#' @param stop_col Stop column name for start-stop Cox model.
#' @param rhs_terms Right-hand side formula terms.
#'
#' @return A formula object.
#' @export

build_surv_formula <- function(event_col,
                               time_col = NULL,
                               start_col = NULL,
                               stop_col = NULL,
                               rhs_terms) {
  rhs_terms <- rhs_terms[!is.na(rhs_terms) & rhs_terms != ""]

  lhs <- if (!is.null(start_col) && !is.null(stop_col)) {
    paste0("survival::Surv(", start_col, ", ", stop_col, ", ", event_col, ")")
  } else {
    paste0("survival::Surv(", time_col, ", ", event_col, ")")
  }

  stats::as.formula(
    paste(lhs, "~", paste(rhs_terms, collapse = " + "))
  )
}


