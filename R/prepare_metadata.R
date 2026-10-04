#' Build baseline or time varying survival metadata from long panel data
#'
#' For baseline survival analysis:
#' Time is set to the first event time if an event occurs, otherwise the last observed follow-up time.
#' For time varying survival analysis:
#' Converts visit-level longitudinal data into counting-process style
#'tart-stop intervals for time-varying Cox analysis.
#'
#' @importFrom rlang .data
#' @param metadata A data frame in long format.
#' @param id_col Subject ID column name.
#' @param time_col Visit time column name.
#' @param event_col Event indicator column name coded 0/1.
#' @param new_time_col Visit time column name created for baseline analysis.
#' @param new_event_col Event indicator column name coded 0/1 created for baseline analysis.
#' @param start_col Start indicator column name
#' @param stop_col Stop indicator column name
#' @param biomarker_type One of "baseline", "time_varying",  or "baseline_change".
#'
#' @details `event_col` must contain non-missing numeric values coded as 0 or 1.
#'   Missing or invalid event statuses cause an error rather than being treated
#'   as censored observations.
#'
#' @details Visit times must be finite, non-missing numeric values with unique subject-time pairs. Longitudinal modes reject first-visit events and require stop > start.
#' @details Output names must be nonempty and distinct, including the subject
#'   ID. Input column names may be reused safely for the output.
#' @return A data frame with one row per subject for baseline mode, or one row
#'   per survival interval for longitudinal modes.


construct_metadata<-function(metadata,
			     id_col = "subject",
                             time_col = "Visit",
                             event_col = NULL,
			     new_time_col='time',
			     new_event_col='event',
			     start_col='start',
			     stop_col='stop',
			     biomarker_type = c("baseline", "time_varying","baseline_change")){

        biomarker_type <- match.arg(biomarker_type)
        if (is.null(event_col) || length(event_col) != 1L ||
            is.na(event_col) || !nzchar(event_col)) {
                stop("`event_col` must name a column in `metadata`.", call. = FALSE)
        }
        if (!is.data.frame(metadata) || !event_col %in% names(metadata)) {
                stop("`event_col` must name a column in `metadata`.", call. = FALSE)
        }
        event_status <- metadata[[event_col]]
        if (!is.numeric(event_status) || anyNA(event_status) ||
            any(!is.finite(event_status)) || any(!event_status %in% c(0, 1))) {
                stop(
                        "`event_col` must contain only non-missing numeric values coded as 0 or 1.",
                        call. = FALSE
                )
        }

        required <- c(id_col, time_col)
        if (!all(required %in% names(metadata))) {
            stop("Subject ID and visit time columns must exist in `metadata`.", call. = FALSE)
        }
        if (anyNA(metadata[[id_col]])) {
            stop("Subject IDs must not be missing.", call. = FALSE)
        }
        times <- metadata[[time_col]]
        if (!is.numeric(times) || anyNA(times) || any(!is.finite(times))) {
            stop("Visit times must be non-missing, finite numeric values.", call. = FALSE)
        }
        output_fields <- if (biomarker_type == "baseline") {
            list(id_col, new_time_col, new_event_col)
        } else list(id_col, start_col, stop_col, event_col)
        if (!all(vapply(output_fields, function(x) is.character(x) &&
                        length(x) == 1L && !is.na(x) && nzchar(x), logical(1)))) {
            stop("Output column names must be nonempty character scalars.", call. = FALSE)
        }
        output_names <- if (biomarker_type == "baseline") {
            c(id_col, new_time_col, new_event_col)
        } else c(id_col, start_col, stop_col, event_col)
        if (anyNA(output_names) || any(!nzchar(output_names)) ||
            anyDuplicated(output_names)) {
            stop("Output column names must be nonempty and distinct, including the subject ID.",
                 call. = FALSE)
        }
        validate_unique_keys(metadata, c(id_col, time_col), "metadata")
        # Normalize a private copy so output names cannot overwrite source fields.
        panel <- metadata[, c(id_col, time_col, event_col), drop = FALSE]
        names(panel) <- c(".subject", ".visit", ".status")
        panel <- dplyr::arrange(panel, .data$.subject, .data$.visit)
        if (biomarker_type != "baseline") {
            first <- !duplicated(panel$.subject)
            affected <- panel$.subject[first & panel$.status == 1]
            if (length(affected)) {
                stop(paste0("First-visit events are unsupported for longitudinal analysis. Subject IDs: ",
                            paste(affected, collapse = ", ")), call. = FALSE)
            }
        }
        if (biomarker_type == "baseline") {
            metadata_format <- panel |>
                dplyr::group_by(.data$.subject) |>
                dplyr::summarise(
                    .survival_time = if (any(.data$.status == 1)) {
                        min(.data$.visit[.data$.status == 1])
                    } else max(.data$.visit),
                    .survival_event = as.integer(any(.data$.status == 1)),
                    .groups = "drop")
        } else {
            metadata_format <- panel |>
                dplyr::group_by(.data$.subject) |>
                dplyr::mutate(.start = .data$.visit,
                              .stop = dplyr::lead(.data$.visit),
                              .interval_event = dplyr::lead(.data$.status)) |>
                dplyr::filter({
                    first_event_row <- match(TRUE, .data$.interval_event == 1)
                    if (is.na(first_event_row)) TRUE else dplyr::row_number() <= first_event_row
                }) |>
                dplyr::filter(!is.na(.data$.stop)) |>
                dplyr::ungroup() |>
                dplyr::select(dplyr::all_of(c(".subject", ".start", ".stop", ".interval_event")))
            metadata_format$.interval_event <- as.integer(metadata_format$.interval_event)
            if (any(metadata_format$.stop <= metadata_format$.start)) {
                stop("Every survival interval must have stop > start.", call. = FALSE)
            }
        }
        names(metadata_format) <- output_names
        metadata_format
}
