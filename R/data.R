#' Toy longitudinal survival metadata
#'
#' Synthetic visit records for four subjects, intended for workflow examples
#' rather than inference. Follow-up time uses arbitrary units.
#' @format A data frame with 10 rows and three columns:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{time}{Visit time in arbitrary follow-up units.}
#'   \item{event}{Event status at the visit: 0 for no event, 1 for an event.}
#' }
"toy_metadata"

#' Toy baseline biomarker measurements
#'
#' Synthetic baseline measurements for five subjects, for workflow examples.
#' Biomarker values use arbitrary measurement units.
#' @format A data frame with five rows and three columns:
#' \describe{
#'   \item{id}{Subject identifier linking to the metadata.}
#'   \item{biomarker1}{First baseline biomarker measurement.}
#'   \item{biomarker2}{Second baseline biomarker measurement.}
#' }
"toy_baseline_biomarkers"

#' Toy longitudinal biomarker measurements
#'
#' Synthetic visit measurements for four subjects, for workflow examples.
#' Time and biomarker values use arbitrary units; time matches toy_metadata.
#' @format A data frame with 10 rows and four columns:
#' \describe{
#'   \item{id}{Subject identifier.}
#'   \item{time}{Visit time in arbitrary follow-up units.}
#'   \item{biomarker1}{First biomarker measurement at the visit.}
#'   \item{biomarker2}{Second biomarker measurement at the visit.}
#' }
"toy_timevarying_biomarkers"

#' Example biomarker selection table
#'
#' Names of the two synthetic biomarkers used in the example workflows.
#' @format A data frame with one row and one column:
#' \describe{
#'   \item{biomarker1}{The value biomarker2. The header names the first
#'     biomarker and the value names the second; use
#'     c("biomarker1", "biomarker2") as the analysis vector.}
#' }
"biomarkers_to_test"
