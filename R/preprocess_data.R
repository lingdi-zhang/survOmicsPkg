#' preproccess survival metadata, biomarker table and covariate table
#'
#' @param metadata Survival metadata,post formatting.
#' @param biomarkers Vector of Biomarker base names.
#' @param id_col Subject ID column name
#' @param time_col Time column name
#' @param new_time_col Visit time column name created for baseline analysis.
#' @param new_event_col Event indicator column name coded 0/1 created for baseline analysis.
#' @param start_col Start column name
#' @param event_col event column name
#' @param stop_col Stop column name
#' @param baseline_feature_table Optional data input for baseline biomarker measurements.
#' @param time_varying_feature_table Optional data input for longitudinal  biomarker measurements.baseline or time varying, need to input either one or both. 
#' @param covariates_table Optional data input for covariates  measurements.
#' @param biomarker_type Biomarker mode.
#' @param baseline_suffix Baseline suffix.
#' @param tv_suffix Time-varying suffix.
#' @param change_suffix Change suffix.
#'
#' @details Rows are retained only when survival fields, selected biomarker
#'   fields, and supplied covariate fields are complete. Missing values in
#'   unused feature columns do not remove rows.
#'
#' @return A data frame with survival metadata and biomarker measurements
#' @export


preprocess_data<-function(metadata,
			   biomarkers,
			   id_col = "subject",
			   time_col="Visit",
			   new_time_col="time",
			   start_col= 'start',
			   event_col=NULL,
			   new_event_col='event',
			   stop_col='stop',
			   baseline_feature_table=NULL,
			   time_varying_feature_table=NULL,
			   covariates_table=NULL,
			   biomarker_type = c("baseline", "time_varying", "baseline_change"),
			   baseline_suffix = "_bl",
			   tv_suffix = "_tv",
			   change_suffix="_delta"){

	biomarker_type <- match.arg(biomarker_type)
	metadata_format=construct_metadata(metadata,
					   id_col = id_col,
                             time_col = time_col,
                             event_col = event_col,
			     new_time_col=new_time_col,
			     new_event_col=new_event_col,
                             start_col=start_col,
                             stop_col=stop_col,
			     biomarker_type=biomarker_type)


	if (biomarker_type =="baseline"){
		baseline_feature_table  <-  baseline_feature_table |>  dplyr::rename_with(~ paste0(.x, baseline_suffix),dplyr::all_of(biomarkers))
		data=dplyr::left_join(metadata_format,baseline_feature_table,by=id_col)
			   }
	if (biomarker_type =="time_varying") {
		time_varying_feature_table  <- time_varying_feature_table |> dplyr::rename_with(~ paste0(.x, tv_suffix),dplyr::all_of(biomarkers))
		time_varying_feature_table <- time_varying_feature_table |>  dplyr::rename(!!start_col := !!time_col)
		keys=c(id_col,start_col)
		data=dplyr::left_join(metadata_format,time_varying_feature_table,by=keys)
	}
	if (biomarker_type =="baseline_change") {
		baseline_feature_table  <- baseline_feature_table |> dplyr::rename_with(~ paste0(.x, baseline_suffix),dplyr::all_of(biomarkers))
#                baseline_feature_table <- baseline_feature_table |>  dplyr::rename(!!start_col := !!time_col)
                time_varying_feature_table  <- time_varying_feature_table |> dplyr::rename_with(~ paste0(.x, tv_suffix),dplyr::all_of(biomarkers))
                time_varying_feature_table <- time_varying_feature_table |>  dplyr::rename(!!start_col := !!time_col)

		keys=c(id_col,start_col)
                data=metadata_format  |>  dplyr::left_join(time_varying_feature_table, by = keys) |>
    		dplyr::left_join(baseline_feature_table, by = id_col)
		for (p in biomarkers) {
            bl_col <- paste0(p, baseline_suffix)
            tv_col <- paste0(p, tv_suffix)
            delta_col <- paste0(p, change_suffix)
            data[[delta_col]] <- data[[tv_col]] - data[[bl_col]]
        }
	
	
	}

	if (!is.null(covariates_table)){
		data=dplyr::inner_join(data,covariates_table,by=id_col)
}

	biomarker_columns <- switch(
		biomarker_type,
		baseline = paste0(biomarkers, baseline_suffix),
		time_varying = paste0(biomarkers, tv_suffix),
		baseline_change = c(
			paste0(biomarkers, baseline_suffix),
			paste0(biomarkers, change_suffix)
		)
	)
	survival_columns <- if (biomarker_type == "baseline") {
		c(id_col, new_time_col, new_event_col)
	} else {
		c(id_col, start_col, stop_col, event_col)
	}
	covariate_columns <- if (is.null(covariates_table)) {
		character()
	} else {
		setdiff(names(covariates_table), id_col)
	}
	required_columns <- unique(c(
		survival_columns,
		biomarker_columns,
		covariate_columns
	))
	data <- data[stats::complete.cases(data[, required_columns, drop = FALSE]), ]
	data
}

	
