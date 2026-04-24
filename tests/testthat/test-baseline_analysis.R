test_that("test baseline analysis", {
	  
   metadata <- data.frame(
    subject = c(1,1,2,2),
    time = c(0, 6, 0, 6),
    event = c(0, 1, 0, 0))
   baseline_df<-data.frame(
	subject = c(1,2),
	biomarker = c(2, 3))

   biomarkers=c("biomarker")
    
  out <- preproccess_data(metadata,biomarkers=biomarkers,event_col="event",time_col="time",baseline_feature_table=baseline_df,biomarker_type="baseline")

  expect_equal(nrow(out), 2)
  expect_true(all(out$event %in% c(0, 1)))

res<-run_multiple_cox_flexible(out,biomarkers=biomarkers,biomarker_type="baseline",time_col="time",event_col="event")
print(res)
outcome=res[[1]]
term_info=res[[2]]
output=calculate_FDR(outcome,term_info)

})




