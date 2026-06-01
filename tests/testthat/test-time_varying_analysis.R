test_that("test time varying analysis", {

   		  
  metadata <- data.frame(
    subject = c(1,1,1,2,2,2),
    time = c(0, 3,6, 0,3, 6),
    event = c(0, 1,1, 0, 0,1))

  time_varying_df<-data.frame(
        subject = c(1,1,1,2,2,2),
        time = c(0, 3,6, 0, 3,6),
        biomarker = c(1,2, 3,4,5,6))

   biomarkers=c("biomarker")


 out <- preprocess_data(metadata,biomarkers=biomarkers,time_col="time",event_col='event',time_varying_feature_table=time_varying_df,biomarker_type="time_varying")

  expect_equal(nrow(out), 3)
  expect_true(all(out$event %in% c(0, 1)))


 res<-run_multiple_cox_flexible(out,biomarker=biomarkers,biomarker_type="time_varying",start_col="start",stop_col="stop",event_col="event")
 print(res)
outcome=res[[1]]
term_info=res[[2]]
output=calculate_FDR(outcome,term_info)

})




