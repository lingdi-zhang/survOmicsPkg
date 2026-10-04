#' Build survival formula
#'
#' @param event_col Event column name.
#' @param time_col Time column name for standard Cox model.
#' @param start_col Start column name for start-stop Cox model.
#' @param stop_col Stop column name for start-stop Cox model.
#' @param rhs_terms Right-hand side formula terms.
#'
#' @details Column names are represented as symbols, supporting spaces and punctuation.
#'   Interval columns must be supplied together.
#' @return A formula object.


build_surv_formula <- function(event_col,
                               time_col = NULL,
                               start_col = NULL,
                               stop_col = NULL,
                               rhs_terms) {
  if (is.null(start_col) != is.null(stop_col)) {
    stop("Supply both `start_col` and `stop_col` together.", call. = FALSE)
  }
  if (is.character(rhs_terms)) rhs_terms <- lapply(rhs_terms, as.name)
  surv_function <- call("::", as.name("survival"), as.name("Surv"))
  lhs <- if (!is.null(start_col) && !is.null(stop_col)) {
    as.call(list(surv_function, as.name(start_col), as.name(stop_col), as.name(event_col)))
  } else {
    if (is.null(time_col)) stop("Supply `time_col` or both interval columns.", call. = FALSE)
    as.call(list(surv_function, as.name(time_col), as.name(event_col)))
  }
  rhs <- if (length(rhs_terms)) Reduce(function(x, y) call("+", x, y), rhs_terms) else 1
  stats::as.formula(call("~", lhs, rhs), env = parent.frame())
}
