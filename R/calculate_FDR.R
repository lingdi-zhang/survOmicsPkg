#' calculate FDRs
#'
#' @param res output from multiple cox models
#' @param terms Optional term_info for legacy results without explicit FDR groups.
#' @param effect_col effect_type column from res.
#' @param term_col Term column from res.
#' @param pvalue_col Numeric p-value column from res; observed values must be
#'   finite and in [0, 1]. Missing values are allowed.
#' @param main_term Main terms from term info
#' @param other_terms Other terms from term info
#' @param other_terms_vector Other_terms or covariate in the effect_type column 
#'
#' @details BH correction is applied within explicit `fdr_group` values.
#'   All rows are preserved; ungrouped, failed, and unreliable fits have missing
#'   FDRs and are excluded from the correction family. Estimable effects from
#'   partial fits remain eligible; not_estimable effects are excluded.
#' @details Legacy grouping matches structural suffixes after unwrapping time coefficients or splitting interactions. Ambiguous matches cause an error; supply explicit fdr_group values.
#' @return A data frame with an FDR column and the original rows preserved.
#' @export

calculate_FDR <- function(res, terms = NULL, effect_col = "effect_type",
                          term_col = "term", pvalue_col = "pvalue",
                          main_term = "main_terms", other_terms = "other_terms",
                          other_terms_vector = c("other_terms", "covariate")) {
  if (!is.data.frame(res) || !is.character(pvalue_col) ||
      length(pvalue_col) != 1L || is.na(pvalue_col) ||
      !pvalue_col %in% names(res)) {
    stop("`pvalue_col` must name an existing column in the result data frame.", call. = FALSE)
  }
  pvalues <- res[[pvalue_col]]
  if (!is.numeric(pvalues)) {
    stop("P-values must be numeric values in [0, 1] or missing.", call. = FALSE)
  }
  observed <- !is.na(pvalues)
  if (any(!is.finite(pvalues[observed])) ||
      any(pvalues[observed] < 0 | pvalues[observed] > 1)) {
    stop("P-values must be finite values in [0, 1] or missing.", call. = FALSE)
  }
  if (!"fdr_group" %in% names(res)) {
    # Legacy results match structural suffixes, assigning each row once.
    if (is.null(terms)) {
      stop("Supply `terms` for results without an `fdr_group` column.", call. = FALSE)
    }
    res$fdr_group <- vapply(seq_len(nrow(res)), function(i) {
      effect <- res[[effect_col]][i]
      term <- res[[term_col]][i]
      if (is.na(effect) || is.na(term)) return(NA_character_)
      is_main <- effect == main_term
      if (!is_main && !effect %in% other_terms_vector) return(NA_character_)
      is_tt <- startsWith(term, "tt(") && endsWith(term, ")")
      structural_term <- if (is_tt) substr(term, 4L, nchar(term) - 1L) else term
      # Parse coefficient labels to respect backticks around nonstandard names.
      expression <- tryCatch(str2lang(structural_term), error = function(e) NULL)
      components <- if (is.symbol(expression)) as.character(expression) else
        if (is.call(expression) && identical(expression[[1L]], as.name(":"))) {
          vapply(as.list(expression)[-1L], as.character, character(1))
        } else character()
      candidates <- terms[[if (is_main || is_tt) main_term else other_terms]]
      matches <- candidates[vapply(candidates, function(candidate) {
        parts <- strsplit(candidate, ":", fixed = TRUE)[[1L]]
        if (length(parts) == 1L) {
          length(components) == 1L && endsWith(components[1L], parts[1L])
        } else if (length(parts) == 2L && length(components) == 2L) {
          (endsWith(components[1L], parts[1L]) && startsWith(components[2L], parts[2L])) ||
            (endsWith(components[2L], parts[1L]) && startsWith(components[1L], parts[2L]))
        } else FALSE
      }, logical(1))]
      if (!length(matches)) return(NA_character_)
      if (length(matches) > 1L) {
        stop(sprintf("Ambiguous FDR group for term `%s`; supply explicit `fdr_group` values.", term),
             call. = FALSE)
      }
      prefix <- if (is_main) "main:" else if (is_tt) "time:" else "other:"
      paste0(prefix, matches)
    }, character(1))
  }
  res$FDR <- rep(NA_real_, nrow(res))
  eligible <- if ("fit_status" %in% names(res)) {
    !is.na(res$fit_status) & res$fit_status %in% c("ok", "warning", "partial")
  } else rep(TRUE, nrow(res))
  if ("effect_status" %in% names(res)) {
    eligible <- eligible & !is.na(res$effect_status) & res$effect_status == "estimable"
  }
  for (group in unique(res$fdr_group[!is.na(res$fdr_group)])) {
    rows <- which(!is.na(res$fdr_group) & res$fdr_group == group & eligible)
    res$FDR[rows] <- stats::p.adjust(res[[pvalue_col]][rows], method = "BH")
  }
  res
}
