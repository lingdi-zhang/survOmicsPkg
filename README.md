# survOmicsPkg



A reproducible framework for biomarker discovery in longitudinal omics studies.

survOmicsPkg is an R package designed to evaluate associations between molecular features and time-to-event outcomes in longitudinal studies. The package implements multiple survival modeling strategies commonly used in translational research, enabling investigators to distinguish prognostic biomarkers, dynamic disease processes, and time-dependent molecular effects within a unified analytical framework.

The framework was developed to support high-dimensional biomarker discovery workflows across transcriptomics, proteomics, metabolomics, microbiome, and clinical datasets.

## Scientific Motivation

Longitudinal biomedical studies frequently collect repeated molecular measurements together with clinical outcomes.

While standard survival analysis tools can fit individual Cox models, applying multiple modeling strategies consistently across large biomarker panels often requires substantial custom code.

survOmicsPkg was developed to provide a reproducible framework for evaluating distinct biological hypotheses, including:

- Whether baseline biomarker levels predict future outcomes
- Whether longitudinal biomarker changes are associated with risk
- Whether biomarker effects vary over time
- Whether risk is driven by baseline differences or within-subject changes

## Analysis Framework

The package supports four complementary survival modeling strategies commonly used in biomarker discovery and translational research:

| Analysis Type | Biological Question |
|--------------|---------------------|
| Baseline Biomarker | Do baseline molecular measurements predict future outcomes? |
| Time-Varying Biomarker | Are longitudinal biomarker changes associated with risk? |
| Time-Varying Coefficient | Does biomarker importance change over time? |
| Baseline + Change | Is risk driven by baseline differences or within-subject changes? |

These complementary modeling approaches allow investigators to evaluate distinct biological hypotheses that are often overlooked in traditional survival analyses. For example, a biomarker may predict risk at baseline but not change over time, while another biomarker may exhibit minimal baseline differences yet become informative through longitudinal changes.

## Workflow
<p align="center">
  <img src="man/figures/Workflow.png" width="850">
</p>

**Figure 1.** Overview of the survOmicsPkg analytical framework. The package supports four complementary survival modeling strategies for evaluating baseline, dynamic, time-dependent, and longitudinal biomarker effects.

## Example Applications

The analytical framework implemented in survOmicsPkg is applicable to:

- Blood transcriptomic studies
- Proteomic analyses
- Metabolomic studies
- Microbiome research
- Longitudinal immune profiling
- Biomarker discovery in translational medicine

The package was designed around common analytical challenges encountered in longitudinal multi-omics studies.


# Installation


``` r
library("remotes")
remotes::install_github("lingdi-zhang/survOmicsPkg", force=TRUE)
library(survOmicsPkg)
```

## Data Requirements

survOmicsPkg requires a longitudinal metadata table and a biomarker dataset linked by a common subject identifier.

### Metadata

| Column | Description |
|----------|-------------|
| id | Subject identifier |
| time | Follow-up time |
| event | Event indicator (0 = censored, 1 = event) |

### Baseline Biomarker Data

| Column | Description |
|----------|-------------|
| id | Subject identifier |
| biomarker1 | Biomarker measurement |
| biomarker2 | Biomarker measurement |
| ... | Additional biomarkers |

### Time-Varying Biomarker Data

| id | time | biomarker1 | biomarker2 |
|----|------|------------|------------|
| 1 | 0 | 5.2 | 3.1 |
| 1 | 2 | 6.8 | 3.4 |
| 1 | 4 | 8.5 | 3.5 |
| 2 | 0 | 4.1 | 2.8 |
| 2 | 2 | 5.9 | 3.0 |
| 2 | 4 | 7.2 | 3.2 |

Each row represents a biomarker measurement collected at a specific follow-up time for a given subject.

### Toy Data for Input Formats

The package includes `toy_metadata`, `toy_baseline_biomarkers`, and
`toy_timevarying_biomarkers`. Use these small tables to inspect input formats
and try preprocessing. The modeling examples below use a larger simulated study.


``` r
data("toy_metadata")
data("toy_baseline_biomarkers")
data("toy_timevarying_biomarkers")
head(toy_metadata)
```

### Simulated Study for Modeling

All four examples use the same 400 subjects and two biomarkers. Follow-up is
measured in years, with annual visits and at most six years of observation.
Event times depend on baseline biomarker values; independent censoring also
occurs. Longitudinal measurements start at the baseline value and vary with
subject-specific slopes and measurement noise. This simulation demonstrates
how to run each model; it does not guarantee significant results for every effect.


``` r
set.seed(2026)
n <- 400
biomarkers <- c("biomarker1", "biomarker2")
sim_baseline <- data.frame(
  id = seq_len(n),
  biomarker1 = rnorm(n),
  biomarker2 = rnorm(n)
)

rate <- 0.15 * exp(
  0.35 * sim_baseline$biomarker1 - 0.25 * sim_baseline$biomarker2
)
event_time <- rexp(n, rate = rate)
censor_time <- rexp(n, rate = 0.03)
followup <- pmin(event_time, censor_time, 6)
event <- as.integer(event_time <= pmin(censor_time, 6))
slopes <- matrix(rnorm(2 * n, sd = 0.2), ncol = 2)

# Include outcome/censoring records as well as all preceding annual visits.
sim_metadata <- do.call(rbind, lapply(seq_len(n), function(i) {
  visits <- c((0:5)[(0:5) < followup[i]], followup[i])
  data.frame(
    id = i,
    time = visits,
    event = c(rep(0L, length(visits) - 1L), event[i])
  )
}))

sim_longitudinal <- sim_metadata[c("id", "time")]
for (j in seq_along(biomarkers)) {
  id <- sim_longitudinal$id
  time <- sim_longitudinal$time
  sim_longitudinal[[biomarkers[j]]] <-
    sim_baseline[[biomarkers[j]]][id] + slopes[id, j] * time +
    rnorm(length(time), sd = 0.15) * sqrt(time)
}

c(subjects = n, events = sum(event), censored = sum(event == 0))
#> subjects   events censored
#>      400      213      187
```

Preprocessing carries each visit's measurement forward to the next visit or
outcome time. Measurements at the outcome time are not used to predict that
same event. The time-varying and baseline-plus-change models below use these
start-stop intervals.

Each result includes `fit_status`, `effect_status`, and `fdr_group`. FDR correction
uses estimable effects in each family, including valid effects from partial fits;
unreliable fits and unestimable effects receive missing FDRs.

## Baseline Biomarker Analysis

Tests whether baseline biomarker levels are associated with future event risk.


``` r
processed_data <- preprocess_data(
  metadata = sim_metadata, biomarkers = biomarkers,
  id_col = "id", event_col = "event", time_col = "time",
  baseline_feature_table = sim_baseline,
  biomarker_type = "baseline"
)
res <- run_multiple_cox_flexible(
  data = processed_data, biomarkers = biomarkers, event_col = "event",
  time_col = "time",
  biomarker_type = "baseline"
)
FDR_res <- calculate_FDR(res[[1]], res[[2]])
FDR_res
#>    biomarker          term       coef        HR   lower95   upper95
#> 1 biomarker1 biomarker1_bl  0.3904288 1.4776142 1.2823335 1.7026334
#> 2 biomarker2 biomarker2_bl -0.2444816 0.7831104 0.6837313 0.8969341
#>         pvalue                                     formula fit_status warning
#> 1 6.719146e-08 survival::Surv(time, event) ~ biomarker1_bl         ok    <NA>
#> 2 4.141333e-04 survival::Surv(time, event) ~ biomarker2_bl         ok    <NA>
#>   effect_type fdr_group effect_status          FDR
#> 1  main_terms  baseline     estimable 1.343829e-07
#> 2  main_terms  baseline     estimable 4.141333e-04
```

**Interpretation**

- HR > 1: Higher baseline biomarker levels are associated with an increased risk of the outcome.
- HR < 1: Higher baseline biomarker levels are associated with a decreased risk of the outcome.

Typical applications include:

- Prognostic biomarker discovery
- Patient risk stratification
- Baseline disease prediction


## Time-Varying Biomarker Analysis

Tests whether biomarker values measured during follow-up are associated with event risk.
This model treats biomarker measurements as time-dependent covariates and updates risk estimates as biomarker values change over follow-up.



``` r
processed_data <- preprocess_data(
  metadata = sim_metadata, biomarkers = biomarkers,
  id_col = "id", event_col = "event", time_col = "time",
  time_varying_feature_table = sim_longitudinal,
  biomarker_type = "time_varying"
)
res <- run_multiple_cox_flexible(
  data = processed_data, biomarkers = biomarkers, event_col = "event",
  start_col = "start", stop_col = "stop",
  biomarker_type = "time_varying"
)
FDR_res <- calculate_FDR(res[[1]], res[[2]])
FDR_res
#>    biomarker          term       coef        HR  lower95  upper95       pvalue
#> 1 biomarker1 biomarker1_tv  0.3192604 1.3761096 1.212782 1.561433 7.319758e-07
#> 2 biomarker2 biomarker2_tv -0.1726188 0.8414583 0.747668 0.947014 4.198225e-03
#>                                              formula fit_status warning
#> 1 survival::Surv(start, stop, event) ~ biomarker1_tv         ok    <NA>
#> 2 survival::Surv(start, stop, event) ~ biomarker2_tv         ok    <NA>
#>   effect_type    fdr_group effect_status          FDR
#> 1  main_terms time_varying     estimable 1.463952e-06
#> 2  main_terms time_varying     estimable 4.198225e-03
```

**Interpretation**

- HR > 1: Higher biomarker levels are associated with an increased risk of the outcome.
- HR < 1: Higher biomarker levels are associated with a decreased risk of the outcome.


This approach incorporates longitudinal biomarker trajectories directly into the survival model.

Typical applications include:

- Dynamic disease monitoring
- Longitudinal immune profiling
- Treatment response biomarkers

## Time-Varying Coefficient Analysis

Tests whether biomarker effects change over time.

This model evaluates interactions between biomarker levels and follow-up time.



``` r
processed_data <- preprocess_data(
  metadata = sim_metadata, biomarkers = biomarkers,
  id_col = "id", event_col = "event", time_col = "time",
  baseline_feature_table = sim_baseline,
  biomarker_type = "baseline"
)
res <- run_multiple_cox_flexible(
  data = processed_data, biomarkers = biomarkers, event_col = "event",
  time_col = "time",
  time_varying_coefficients = TRUE, center_time = 0,
  biomarker_type = "baseline"
)
FDR_res <- calculate_FDR(res[[1]], res[[2]])
FDR_res
#>    biomarker              term        coef        HR   lower95   upper95
#> 1 biomarker1     biomarker1_bl  0.44144614 1.5549543 1.2214894 1.9794545
#> 2 biomarker1 tt(biomarker1_bl) -0.02217183 0.9780722 0.8985792 1.0645974
#> 3 biomarker2     biomarker2_bl -0.26961492 0.7636735 0.6057325 0.9627966
#> 4 biomarker2 tt(biomarker2_bl)  0.01081805 1.0108768 0.9324568 1.0958919
#>         pvalue                                                         formula
#> 1 0.0003376753 survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 2 0.6082011857 survival::Surv(time, event) ~ biomarker1_bl + tt(biomarker1_bl)
#> 3 0.0225679573 survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#> 4 0.7928789384 survival::Surv(time, event) ~ biomarker2_bl + tt(biomarker2_bl)
#>   fit_status warning effect_type        fdr_group effect_status          FDR
#> 1         ok    <NA>  main_terms         baseline     estimable 0.0006753505
#> 2         ok    <NA> other_terms time_coefficient     estimable 0.7928789384
#> 3         ok    <NA>  main_terms         baseline     estimable 0.0225679573
#> 4         ok    <NA> other_terms time_coefficient     estimable 0.7928789384
```

**Interpretation**

- biomarker*_bl:
Tests the baseline biomarker association at `center_time` (year zero here).
- tt(biomarker*_bl):
Tests whether the association between baseline biomarker levels and the outcome changes over time.

For biomarker*_bl, evaluated at `center_time`:

- HR > 1: Higher baseline biomarker levels are associated with an increased risk of the outcome.
- HR < 1: Higher baseline biomarker levels are associated with a decreased risk of the outcome.

The hazard ratio at year `t` is `exp(beta_bl + beta_tt * (t - center_time))`.
The exponentiated `tt()` coefficient describes the change in this hazard ratio
per year, rather than a separate overall biomarker hazard ratio.

For tt(biomarker*_bl):

- HR > 1: The effect of the baseline biomarker on outcome risk increases over time.
- HR < 1: The effect of the baseline biomarker on outcome risk decreases over time.

**Note:** The simulated event process uses constant baseline effects. A time-varying coefficient can therefore be estimable without providing evidence that the effect changes over time.

Typical applications include:

- Early versus late disease effects
- Treatment adaptation studies
- Temporal changes in biomarker importance

## Baseline + Change Model

Separates:

- Between-subject effects (baseline differences)
- Within-subject effects (longitudinal changes)

This decomposition helps determine whether event risk is driven by:

1. Persistent baseline differences between individuals
2. Dynamic biomarker changes within individuals over time

This framework is particularly useful for repeated-measurement studies.



``` r
processed_data <- preprocess_data(
  metadata = sim_metadata, biomarkers = biomarkers,
  id_col = "id", event_col = "event", time_col = "time",
  baseline_feature_table = sim_baseline,
  time_varying_feature_table = sim_longitudinal,
  biomarker_type = "baseline_change"
)
res <- run_multiple_cox_flexible(
  data = processed_data, biomarkers = biomarkers, event_col = "event",
  start_col = "start", stop_col = "stop",
  biomarker_type = "baseline_change"
)
FDR_res <- calculate_FDR(res[[1]], res[[2]])
FDR_res
#>    biomarker             term        coef        HR   lower95   upper95
#> 1 biomarker1    biomarker1_bl  0.39038541 1.4775502 1.2822866 1.7025480
#> 2 biomarker1 biomarker1_delta  0.04615852 1.0472404 0.7933065 1.3824574
#> 3 biomarker2    biomarker2_bl -0.24571286 0.7821468 0.6829794 0.8957131
#> 4 biomarker2 biomarker2_delta  0.08845728 1.0924876 0.8432959 1.4153148
#>         pvalue
#> 1 6.731936e-08
#> 2 7.445947e-01
#> 3 3.821392e-04
#> 4 5.030705e-01
#>                                                                 formula
#> 1 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 2 survival::Surv(start, stop, event) ~ biomarker1_bl + biomarker1_delta
#> 3 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#> 4 survival::Surv(start, stop, event) ~ biomarker2_bl + biomarker2_delta
#>   fit_status warning effect_type fdr_group effect_status          FDR
#> 1         ok    <NA>  main_terms  baseline     estimable 1.346387e-07
#> 2         ok    <NA>  main_terms    change     estimable 7.445947e-01
#> 3         ok    <NA>  main_terms  baseline     estimable 3.821392e-04
#> 4         ok    <NA>  main_terms    change     estimable 7.445947e-01
```
**Interpretation**

- biomarker*_bl:
Tests the association between baseline biomarker levels and the outcome.

- biomarker*_delta:
Tests the association between within-subject changes in biomarker levels and the outcome.

For biomarker*_bl:

- HR > 1: Higher baseline biomarker levels are associated with an increased risk of the outcome.
- HR < 1: Higher baseline biomarker levels are associated with a decreased risk of the outcome.

For biomarker*_delta:

- HR > 1: Larger increases (or smaller decreases) in biomarker levels are associated with an increased risk of the outcome.
- HR < 1: Larger increases (or smaller decreases) in biomarker levels are associated with a decreased risk of the outcome.

## Technical Design

The package was designed to provide a consistent interface for evaluating multiple biological hypotheses within a unified survival modeling framework.

Key design principles include:

- Modular implementation of multiple survival modeling strategies
- Consistent handling of longitudinal and time-to-event data
- Automated extraction of hazard ratios, confidence intervals, and significance statistics
- Support for high-throughput biomarker screening
- Reproducible workflows for large-scale omics analyses

## Repository Structure

```text

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

survOmicsPkg was developed to address analytical challenges commonly encountered in longitudinal multi-omics studies. By providing a unified framework for evaluating baseline, dynamic, and time-dependent biomarker effects, the package facilitates reproducible biomarker discovery and translational research workflows.

## Methods Implemented

The package currently supports:

- Cox proportional hazards models
- Time-varying covariate Cox models
- Time-varying coefficient models
- Baseline-versus-change decomposition
- Multiple-testing correction using FDR

Future development will focus on joint longitudinal-survival models and additional time-to-event methodologies.
