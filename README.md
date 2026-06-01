survOmicsPkg
================

A reproducible framework for biomarker discovery in longitudinal omics
studies.

survOmicsPkg is an R package designed to evaluate associations between
molecular features and time-to-event outcomes in longitudinal studies.
The package implements multiple survival modeling strategies commonly
used in translational research, enabling investigators to distinguish
prognostic biomarkers, dynamic disease processes, and time-dependent
molecular effects within a unified analytical framework.

The framework was developed to support high-dimensional biomarker
discovery workflows across transcriptomics, proteomics, metabolomics,
microbiome, and clinical datasets.

## Scientific Motivation

Longitudinal biomedical studies frequently collect repeated molecular
measurements together with clinical outcomes.

While standard survival analysis tools can fit individual Cox models,
applying multiple modeling strategies consistently across large
biomarker panels often requires substantial custom code.

survOmicsPkg was developed to provide a reproducible framework for
evaluating distinct biological hypotheses, including:

-   Whether baseline biomarker levels predict future outcomes
-   Whether longitudinal biomarker changes are associated with risk
-   Whether biomarker effects vary over time
-   Whether risk is driven by baseline differences or within-subject
    changes

## Analysis Framework

The package supports four complementary survival modeling strategies
commonly used in biomarker discovery and translational research:

| Analysis Type            | Biological Question                                               |
|--------------------------|-------------------------------------------------------------------|
| Baseline Biomarker       | Do baseline molecular measurements predict future outcomes?       |
| Time-Varying Biomarker   | Are longitudinal biomarker changes associated with risk?          |
| Time-Varying Coefficient | Does biomarker importance change over time?                       |
| Baseline + Change        | Is risk driven by baseline differences or within-subject changes? |

These complementary modeling approaches allow investigators to evaluate
distinct biological hypotheses that are often overlooked in traditional
survival analyses. For example, a biomarker may predict risk at baseline
but not change over time, while another biomarker may exhibit minimal
baseline differences yet become informative through longitudinal
changes.

## Workflow

<p align="center">
<img src="man/figures/Workflow.png" width="850">
</p>

**Figure 1.** Overview of the survOmicsPkg analytical framework. The
package supports four complementary survival modeling strategies for
evaluating baseline, dynamic, time-dependent, and longitudinal biomarker
effects.

## Example Applications

The analytical framework implemented in survOmicsPkg is applicable to:

-   Blood transcriptomic studies
-   Proteomic analyses
-   Metabolomic studies
-   Microbiome research
-   Longitudinal immune profiling
-   Biomarker discovery in translational medicine

The package was designed around common analytical challenges encountered
in longitudinal multi-omics studies.

# Installation

``` r
library("remotes")
remotes::install_github("lingdi-zhang/survOmicsPkg", force=TRUE)
#> Using GitHub PAT from the git credential store.
#> Downloading GitHub repo lingdi-zhang/survOmicsPkg@HEAD
#> 
#> ── R CMD build ─────────────────────────────────────────────────────────────────
#> * checking for file ‘/private/var/folders/0x/ctzdrfzd3p3b27z358ztmgyc0000gn/T/RtmpKyAe8j/remotes6779495dc694/lingdi-zhang-survOmicsPkg-0f264d4/DESCRIPTION’ ... OK
#> * preparing ‘survOmicsPkg’:
#> * checking DESCRIPTION meta-information ... OK
#> * excluding invalid files
#> Subdirectory 'man' contains invalid file names:
#>   ‘Workflow.png’
#> * checking for LF line-endings in source and make files and shell scripts
#> * checking for empty or unneeded directories
#> * building ‘survOmicsPkg_0.1.0.tar.gz’
library(survOmicsPkg)
```

## Data Requirements

survOmicsPkg requires a longitudinal metadata table and a biomarker
dataset linked by a common subject identifier.

### Metadata

| Column | Description                               |
|--------|-------------------------------------------|
| id     | Subject identifier                        |
| time   | Follow-up time                            |
| event  | Event indicator (0 = censored, 1 = event) |

### Baseline Biomarker Data

| Column     | Description           |
|------------|-----------------------|
| id         | Subject identifier    |
| biomarker1 | Biomarker measurement |
| biomarker2 | Biomarker measurement |
| …          | Additional biomarkers |

### Time-Varying Biomarker Data

| id  | time | biomarker1 | biomarker2 |
|-----|------|------------|------------|
| 1   | 0    | 5.2        | 3.1        |
| 1   | 2    | 6.8        | 3.4        |
| 1   | 4    | 8.5        | 3.5        |
| 2   | 0    | 4.1        | 2.8        |
| 2   | 2    | 5.9        | 3.0        |
| 2   | 4    | 7.2        | 3.2        |

Each row represents a biomarker measurement collected at a specific
follow-up time for a given subject.

### Example Datasets

``` r
data("toy_metadata")
data("toy_baseline_biomarkers")
data("toy_timevarying_biomarkers")
biomarkers <- c("biomarker1", "biomarker2")
```

## Baseline Biomarker Analysis

Tests whether baseline biomarker levels are associated with future event
risk.

``` r
processed_data<-preprocess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",baseline_feature_table=toy_baseline_biomarkers,biomarker_type="baseline")

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",time_col="time",biomarker_type="baseline")

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

FDR_res
#>    biomarker          term       coef        HR   lower95  upper95    pvalue
#> 1 biomarker1 biomarker1_bl  0.2110877 1.2350207 0.6096223 2.502002 0.5578696
#> 2 biomarker2 biomarker2_bl -0.3706716 0.6902706 0.2268980 2.099946 0.5137642
#>                                       formula effect_type       FDR
#> 1 survival::Surv(time, event) ~ biomarker1_bl  main_terms 0.5578696
#> 2 survival::Surv(time, event) ~ biomarker2_bl  main_terms 0.5578696
```

**Interpretation**

-   HR &gt; 1: Higher baseline biomarker levels are associated with an
    increased risk of the outcome.
-   HR &lt; 1: Higher baseline biomarker levels are associated with a
    decreased risk of the outcome.

Typical applications include:

-   Prognostic biomarker discovery
-   Patient risk stratification
-   Baseline disease prediction

## Time-Varying Biomarker Analysis

Tests whether biomarker values measured during follow-up are associated
with event risk. This model treats biomarker measurements as
time-dependent covariates and updates risk estimates as biomarker values
change over follow-up.

``` r
processed_data<-preprocess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",start_col= 'start', stop_col='stop',time_varying_feature_table=toy_timevarying_biomarkers,biomarker_type="time_varying")

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",start_col= 'start', stop_col='stop',biomarker_type="time_varying")

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

FDR_res
#>    biomarker          term       coef        HR   lower95  upper95    pvalue
#> 1 biomarker1 biomarker1_tv  0.5676601 1.7641344 0.4039194 7.704929 0.4504231
#> 2 biomarker2 biomarker2_tv -0.7563076 0.4693964 0.0678334 3.248149 0.4434946
#>                                              formula effect_type       FDR
#> 1 survival::Surv(start, stop, event) ~ biomarker1_tv  main_terms 0.4504231
#> 2 survival::Surv(start, stop, event) ~ biomarker2_tv  main_terms 0.4504231
```

**Interpretation**

-   HR &gt; 1: Higher biomarker levels are associated with an increased
    risk of the outcome.
-   HR &lt; 1: Higher biomarker levels are associated with a decreased
    risk of the outcome.

This approach incorporates longitudinal biomarker trajectories directly
into the survival model.

Typical applications include:

-   Dynamic disease monitoring
-   Longitudinal immune profiling
-   Treatment response biomarkers

## Time-Varying Coefficient Analysis

Tests whether biomarker effects change over time.

This model evaluates interactions between biomarker levels and follow-up
time.

``` r
processed_data<-preprocess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",baseline_feature_table=toy_baseline_biomarkers,biomarker_type="baseline")

res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",time_col="time",biomarker_type="baseline",time_varying_coefficients=TRUE,center_time=0)

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

FDR_res
#>       biomarker              term       coef        HR   lower95  upper95
#> 1 biomarker1_bl     biomarker1_bl  0.2110877 1.2350207 0.6096223 2.502002
#> 2 biomarker2_bl     biomarker2_bl -0.3706716 0.6902706 0.2268980 2.099946
#> 3 biomarker1_bl tt(biomarker1_bl)         NA        NA        NA       NA
#> 4 biomarker2_bl tt(biomarker2_bl)         NA        NA        NA       NA
#>      pvalue                                                         formula
#> 1 0.5578696 survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 2 0.5137642 survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#> 3        NA survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 4        NA survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#>   effect_type       FDR
#> 1  main_terms 0.5578696
#> 2  main_terms 0.5578696
#> 3   covariate        NA
#> 4   covariate        NA
```

**Interpretation**

-   biomarker\*\_bl: Tests the association between baseline biomarker
    levels and the outcome.
-   tt(biomarker\*\_bl): Tests whether the association between baseline
    biomarker levels and the outcome changes over time.

For biomarker\*\_bl:

-   HR &gt; 1: Higher baseline biomarker levels are associated with an
    increased risk of the outcome.
-   HR &lt; 1: Higher baseline biomarker levels are associated with a
    decreased risk of the outcome.

For tt(biomarker\*\_bl):

-   HR &gt; 1: The effect of the baseline biomarker on outcome risk
    increases over time.
-   HR &lt; 1: The effect of the baseline biomarker on outcome risk
    decreases over time.

**Note:** Time-varying coefficient models typically require larger
sample sizes and more events than standard Cox models. The toy dataset
is provided for demonstration purposes only, and some interaction terms
may not be estimable.

Typical applications include:

-   Early versus late disease effects
-   Treatment adaptation studies
-   Temporal changes in biomarker importance

## Baseline + Change Model

Separates:

-   Between-subject effects (baseline differences)
-   Within-subject effects (longitudinal changes)

This decomposition helps determine whether event risk is driven by:

1.  Persistent baseline differences between individuals
2.  Dynamic biomarker changes within individuals over time

This framework is particularly useful for repeated-measurement studies.

``` r
#Both the baseline biomarkers and time varying biomarkers are required in this model
processed_data<-preprocess_data(
  metadata=toy_metadata,biomarkers=biomarkers,id_col="id",event_col="event",time_col="time",start_col= 'start', stop_col='stop',baseline_feature_table=toy_baseline_biomarkers,time_varying_feature_table=toy_timevarying_biomarkers,biomarker_type="baseline_change")


res<-run_multiple_cox_flexible(processed_data,biomarkers=biomarkers, event_col="event",start_col= 'start', stop_col='stop',biomarker_type="baseline_change")

outcomes=res[[1]]
terms=res[[2]]
FDR_res<-calculate_FDR(outcomes,terms)

FDR_res
#>    biomarker             term        coef        HR    lower95   upper95
#> 1 biomarker1    biomarker1_bl  0.51812371 1.6788746 0.37749434  7.466655
#> 2 biomarker2    biomarker2_bl -0.50376179 0.6042533 0.06744550  5.413587
#> 3 biomarker1 biomarker1_delta  0.30259985 1.3533728 0.15544544 11.783028
#> 4 biomarker2 biomarker2_delta -0.04383398 0.9571128 0.08722177 10.502711
#>      pvalue
#> 1 0.4961974
#> 2 0.6524961
#> 3 0.7840377
#> 4 0.9713901
#>                                                                 formula
#> 1 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 2 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#> 3 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 4 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#>   effect_type       FDR
#> 1  main_terms 0.6524961
#> 2  main_terms 0.6524961
#> 3  main_terms 0.9713901
#> 4  main_terms 0.9713901
```

**Interpretation**

-   biomarker\*\_bl: Tests the association between baseline biomarker
    levels and the outcome.

-   biomarker\*\_delta: Tests the association between within-subject
    changes in biomarker levels and the outcome.

For biomarker\*\_bl:

-   HR &gt; 1: Higher baseline biomarker levels are associated with an
    increased risk of the outcome.
-   HR &lt; 1: Higher baseline biomarker levels are associated with a
    decreased risk of the outcome.

For biomarker\*\_delta:

-   HR &gt; 1: Larger increases (or smaller decreases) in biomarker
    levels are associated with an increased risk of the outcome.
-   HR &lt; 1: Larger increases (or smaller decreases) in biomarker
    levels are associated with a decreased risk of the outcome.

## Technical Design

The package was designed to provide a consistent interface for
evaluating multiple biological hypotheses within a unified survival
modeling framework.

Key design principles include:

-   Modular implementation of multiple survival modeling strategies
-   Consistent handling of longitudinal and time-to-event data
-   Automated extraction of hazard ratios, confidence intervals, and
    significance statistics
-   Support for high-throughput biomarker screening
-   Reproducible workflows for large-scale omics analyses

## Repository Structure

``` text
survOmicsPkg/
├── R/                    # Core package functions
├── man/                  # Function documentation
├── tests/                # Unit tests
├── data/                 # Example datasets
├── man/figures/          # Workflow and documentation figures
├── README.Rmd
├── README.md
├── DESCRIPTION
└── NAMESPACE
```

## Applications in Translational Research

survOmicsPkg was developed to address analytical challenges commonly
encountered in longitudinal multi-omics studies. By providing a unified
framework for evaluating baseline, dynamic, and time-dependent biomarker
effects, the package facilitates reproducible biomarker discovery and
translational research workflows.

## Methods Implemented

The package currently supports:

-   Cox proportional hazards models
-   Time-varying covariate Cox models
-   Time-varying coefficient models
-   Baseline-versus-change decomposition
-   Multiple-testing correction using FDR

Future development will focus on joint longitudinal-survival models and
additional time-to-event methodologies.
