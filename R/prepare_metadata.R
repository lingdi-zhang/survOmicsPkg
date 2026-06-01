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
#' @return A data frame with one row per subject.


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
        if (biomarker_type=="baseline") {
                metadata_format <-
                dplyr::arrange(metadata,.data[[id_col]], .data[[time_col]]) |>
                dplyr::group_by(.data[[id_col]]) |>
  		dplyr::summarise(!!new_event_col:= as.integer(any(.data[[event_col]] == 1)),
                                 !!new_time_col:= ifelse(any(.data[[event_col]] == 1),
                                 min(.data[[time_col]][.data[[event_col]] == 1]),
                                 max(.data[[time_col]])),
                                 .groups = "drop")
	}

     	event_interval <- NULL	
        if (biomarker_type %in% c("time_varying","baseline_change")) {
                metadata_format <- 
                dplyr::arrange(metadata,.data[[id_col]], .data[[time_col]]) |>
                dplyr::group_by(.data[[id_col]]) |>
                dplyr::mutate(!!start_col := .data[[time_col]],
                        !!stop_col := dplyr::lead(.data[[time_col]]),
                        event_interval = dplyr::lead(.data[[event_col]])) |>
		dplyr::filter({first_event_row <- match(TRUE, event_interval == 1)
                       if (is.na(first_event_row)) TRUE else dplyr::row_number() <= first_event_row}) |>
                dplyr::filter(!is.na(.data[[stop_col]])) |>
                dplyr::mutate(!!event_col := as.integer(dplyr::coalesce(event_interval, 0L))) |>
		dplyr::ungroup() |>
                dplyr::select(dplyr::all_of(id_col), dplyr::all_of(start_col), dplyr::all_of(stop_col), dplyr::all_of(event_col))}

	metadata_format
}





