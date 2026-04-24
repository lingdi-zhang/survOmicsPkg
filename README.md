survOmicsPkg
================

survOmicsPkg provides tools for preparing longitudinal data for survival
analysis. It supports construction of metadata for baseline biomarker
and time varying biomarker analyses and flexible survival analyses on
omics datasets.

Install the package from Github:

``` r
library("remotes")
remotes::install_github("lingdi-zhang/survOmicsPkg", force=TRUE)
#> Downloading GitHub repo lingdi-zhang/survOmicsPkg@HEAD
#> Running `R CMD build`...
#> * checking for file ‘/private/var/folders/0x/ctzdrfzd3p3b27z358ztmgyc0000gn/T/RtmpOSnCAL/remotes3f3a54418e4c/lingdi-zhang-survOmicsPkg-19f609e/DESCRIPTION’ ... OK
#> * preparing ‘survOmicsPkg’:
#> * checking DESCRIPTION meta-information ... OK
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘survOmicsPkg_0.0.0.9000.tar.gz’
```

Overview

This package is designed for workflows involving:

longitudinal or repeated-measures data  
survival or time-to-event analysis  
reproducible preprocessing of omics or clinical datasets

A typical workflow:

Prepare longitudinal data in long format  
Construct datasets using preprocess\_data()  
run run\_multiple\_cox\_flexible() for omics data with multiple
biomarkers  
run calculate\_FDR for compute FDR for each effect type (main effect,
time effect, covariates or interactions term) from the results

Load required packages:

``` r
library(survOmicsPkg)
```

Input Metadata Format

The package expects a long-format metadata dataset where each row
represents one subject at one time point.

``` r
toy_metadata <- data.frame(
  id = c(1,1,1,2,2,3,3),
  time = c(0,2,6,0,2,0,2),
  event = c(0,0,1,0,1,0,0))
toy_metadata
#>   id time event
#> 1  1    0     0
#> 2  1    2     0
#> 3  1    6     1
#> 4  2    0     0
#> 5  2    2     1
#> 6  3    0     0
#> 7  3    2     0
```

Required columns for metadata

id: subject identifier  
time: time variable (numeric)  
event: event indicator (0/1)

For baseline biomarker analysis:

The package expects a baseline biomarker dataframe and a vector of
biomarker feature names

``` r
baseline_biomarkers<-data.frame(
  id=c(1,2,3),
  biomarker1=c(10.2,8.5,7.2),
  biomarker2=c(5.1,4.8,6.5)
)
baseline_biomarkers
#>   id biomarker1 biomarker2
#> 1  1       10.2        5.1
#> 2  2        8.5        4.8
#> 3  3        7.2        6.5

biomarkers=c("biomarker1","biomarker2")
biomarkers
#> [1] "biomarker1" "biomarker2"
```

Required columns for baseline biomarker dataset

id: subject identifier  
feature columns: biomarkers, genes, or other measurements

This analysis test for: h(t∣x)=h0(t)exp(βx)  
β is biomarker effect on log hazard ratio

``` r
processed_data<-preproccess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",baseline_feature_table=baseline_biomarkers,biomarker_type="baseline")

processed_data
#> # A tibble: 3 × 5
#>      id event  time biomarker1_bl biomarker2_bl
#>   <dbl> <int> <dbl>         <dbl>         <dbl>
#> 1     1     1     0          10.2           5.1
#> 2     2     1     0           8.5           4.8
#> 3     3     0     2           7.2           6.5

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",time_col="time",biomarker_type="baseline")

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

head(FDR_res)
#>    biomarker          term       coef        HR    lower95   upper95
#> 1 biomarker1 biomarker1_bl  0.5906921 1.8052375 0.53885644  6.047775
#> 2 biomarker2 biomarker2_bl -2.0649225 0.1268281 0.00123326 13.042973
#>      pvalue                                     formula effect_type
#> 1 0.3382651 survival::Surv(time, event) ~ biomarker1_bl  main_terms
#> 2 0.3823791 survival::Surv(time, event) ~ biomarker2_bl  main_terms
#>         FDR
#> 1 0.3823791
#> 2 0.3823791
```

For baseline analysis with time varying outcomes:

set time\_varying\_coefficients=TRUE and adjust center\_time if needed

This analysis test for: h(t∣x)=h0(t)exp(β1x+β2x⋅t)  
β1 is biomarker effect on log hazard ratio at t=0  
β2 is the interaction effect: how the main effect change over time

This is an example, with the toy data, the interaction term won’t
converge.

``` r
res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",time_col="time",biomarker_type="baseline",time_varying_coefficients=TRUE,center_time=0)

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

head(FDR_res)
#>       biomarker              term       coef        HR    lower95
#> 1 biomarker1_bl     biomarker1_bl  0.5906921 1.8052375 0.53885644
#> 2 biomarker2_bl     biomarker2_bl -2.0649225 0.1268281 0.00123326
#> 3 biomarker1_bl tt(biomarker1_bl)         NA        NA         NA
#> 4 biomarker2_bl tt(biomarker2_bl)         NA        NA         NA
#>     upper95    pvalue
#> 1  6.047775 0.3382651
#> 2 13.042973 0.3823791
#> 3        NA        NA
#> 4        NA        NA
#>                                                           formula
#> 1 survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 2 survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#> 3 survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 4 survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#>   effect_type       FDR
#> 1  main_terms 0.3823791
#> 2  main_terms 0.3823791
#> 3   covariate        NA
#> 4   covariate        NA
```

For time varying biomarker analysis: The package expects a time varying
biomarker dataframe and a vector of biomarker feature names

This analysis tests for: h(t∣x(t))=h0(t)exp(β⋅x(t))  
β is biomarker effect on log hazard ratio, while the biomarker values
changes over time

``` r
time_varying_biomarkers<-data.frame(
  id = c(1,1,1,2,2,3,3),
  time = c(0,2,6,0,2,0,2),
  biomarker1 = c(10.2, 11.4, 13.1, 8.5, 9.7, 7.2, 7.8),
  biomarker2 = c(5.1, 5.5, 6.2, 4.8, 5.0, 6.5, 6.3)
)

time_varying_biomarkers
#>   id time biomarker1 biomarker2
#> 1  1    0       10.2        5.1
#> 2  1    2       11.4        5.5
#> 3  1    6       13.1        6.2
#> 4  2    0        8.5        4.8
#> 5  2    2        9.7        5.0
#> 6  3    0        7.2        6.5
#> 7  3    2        7.8        6.3
biomarkers=c("biomarker1","biomarker2")
biomarkers
#> [1] "biomarker1" "biomarker2"
```

Required columns for time varying biomarker dataset

id: subject identifier  
time: time variable (numeric)  
feature columns: biomarkers, genes, or other measurements

``` r
processed_data<-preproccess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",start_col= 'start', stop_col='stop',time_varying_feature_table=time_varying_biomarkers,biomarker_type="time_varying")

processed_data
#> # A tibble: 4 × 6
#>      id start  stop event biomarker1_tv biomarker2_tv
#>   <dbl> <dbl> <dbl> <int>         <dbl>         <dbl>
#> 1     1     0     2     0          10.2           5.1
#> 2     1     2     6     1          11.4           5.5
#> 3     2     0     2     1           8.5           4.8
#> 4     3     0     2     0           7.2           6.5

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",start_col= 'start', stop_col='stop',biomarker_type="time_varying")
#> Warning in agreg.fit(X, Y, istrat, offset, init, control, weights =
#> weights, : Ran out of iterations and did not converge

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

head(FDR_res)
#>    biomarker          term         coef           HR   lower95
#> 1 biomarker1 biomarker1_tv  -0.08942133 9.144602e-01 0.1814907
#> 2 biomarker2 biomarker2_tv -62.68041999 6.001219e-28 0.0000000
#>    upper95    pvalue
#> 1 4.607605 0.9136952
#> 2      Inf 0.9987614
#>                                              formula effect_type
#> 1 survival::Surv(start, stop, event) ~ biomarker1_tv  main_terms
#> 2 survival::Surv(start, stop, event) ~ biomarker2_tv  main_terms
#>         FDR
#> 1 0.9987614
#> 2 0.9987614
```

For analysis consider both baseline biomarker effect and within-subject
change from baseline:

This analysis tests for: h(t)=h0(t)exp(β1x0+β2Δx(t)) β1 represents the
baseline biomarker effect on log hazard  
β2 represents the effect of within-subject change from baseline on log
hazard

This is an exmaple, with the toy dataset, the β2 term won’t converge.

``` r
processed_data<-preproccess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",start_col= 'start', stop_col='stop',baseline_feature_table=baseline_biomarkers,time_varying_feature_table=time_varying_biomarkers,biomarker_type="baseline_change")

processed_data
#> # A tibble: 4 × 10
#>      id start  stop event biomarker1_tv biomarker2_tv biomarker1_bl
#>   <dbl> <dbl> <dbl> <int>         <dbl>         <dbl>         <dbl>
#> 1     1     0     2     0          10.2           5.1          10.2
#> 2     1     2     6     1          11.4           5.5          10.2
#> 3     2     0     2     1           8.5           4.8           8.5
#> 4     3     0     2     0           7.2           6.5           7.2
#> # ℹ 3 more variables: biomarker2_bl <dbl>, biomarker1_delta <dbl>,
#> #   biomarker2_delta <dbl>

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",start_col= 'start', stop_col='stop',biomarker_type="baseline_change")
#> Warning in agreg.fit(X, Y, istrat, offset, init, control, weights =
#> weights, : Ran out of iterations and did not converge

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

head(FDR_res)
#>    biomarker             term         coef           HR   lower95
#> 1 biomarker1    biomarker1_bl  -0.08942133 9.144602e-01 0.1814907
#> 2 biomarker2    biomarker2_bl -52.68039476 1.321890e-23 0.0000000
#> 3 biomarker1 biomarker1_delta           NA           NA        NA
#> 4 biomarker2 biomarker2_delta           NA           NA        NA
#>    upper95    pvalue
#> 1 4.607605 0.9136952
#> 2      Inf 0.9953346
#> 3       NA        NA
#> 4       NA        NA
#>                                                                 formula
#> 1 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 2 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#> 3 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 4 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#>   effect_type       FDR
#> 1  main_terms 0.9953346
#> 2  main_terms 0.9953346
#> 3  main_terms        NA
#> 4  main_terms        NA
```
