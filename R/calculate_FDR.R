#' calculate FDRs
#'
#' @param res output from multiple cox models
#' @param terms term_info from term builders
#' @param effect_col effect_type column from res.
#' @param term_col Term column from res.
#' @param pvalue_col P value column from res
#' @param main_term Main terms from term info
#' @param other_terms Other terms from term info
#' @param other_terms_vector Other_terms or covariate in the effect_type column 
#'
#' @return A data frame with FDR calculated per term
#' @export

calculate_FDR<-function(res,terms,effect_col='effect_type',term_col='term',pvalue_col='pvalue',main_term='main_terms',other_terms="other_terms",other_terms_vector=c('other_terms',"covariate")){
	main_data <- res |> dplyr::filter(.data[[effect_col]] == main_term)
	interaction_data=res[which(res[[effect_col]] %in% other_terms_vector),]
	main_output=list()
	interaction_output=list()
	for (i in (1:length(terms[[main_term]]))){
		main <- main_data |> dplyr::filter(stringr::str_detect(.data[[term_col]], terms[[main_term]][i]))
		main$FDR=stats::p.adjust(main[[pvalue_col]],method='fdr')
		main_output[[i]]=main
	}

	if (length(terms[[other_terms]]) >0) {
	for (i in (1:length(terms[[other_terms]]))){
                interaction <- interaction_data |> dplyr::filter(stringr::str_detect(.data[[term_col]], terms[[other_terms]][i]))
                interaction$FDR=stats::p.adjust(interaction[[pvalue_col]],method='fdr')
                interaction_output[[i]]=interaction
        }

	main_OUT=data.frame(data.table::rbindlist(main_output))
	interaction_OUT=data.frame(data.table::rbindlist(interaction_output))

	OUT=rbind(main_OUT,interaction_OUT)}
	else {
		OUT=data.frame(data.table::rbindlist(main_output))}

	OUT
}






