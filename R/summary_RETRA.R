#' Summarize representative trajectories
#'
#' @description
#' Summarize the properties of representative trajectories returned by
#' [`retra_edr()`] or [`define_retra()`]
#'
#' @param object An object of class `RETRA`.
#' @param ... (not used)
#'
#' @return
#' Data frame with nine columns and one row for each representative trajectory
#' in `object`. The columns in the returned data frame contain the following
#' information:
#' \describe{
#' \item{`ID`}{Identifier of the representative trajectories.}
#' \item{`Size`}{Number of states forming each representative trajectory.}
#' \item{`Length`}{Sum of the dissimilarities in `d` between every pair of
#' consecutive states forming the representative trajectories.}
#' \item{`Avg_link`}{Mean value of the dissimilarities between consecutive states
#' of the representative trajectories that do not belong to the same ecological
#' trajectory or site (i.e., artificial links).}
#' \item{`Sum_link`}{Sum of the dissimilarities between consecutive states of the
#' representative trajectories that do not belong to the same ecological trajectory
#' or site (i.e., artificial links).}
#' \item{`Avg_density`}{Mean value of the number of segments represented by each
#' segment of the representative trajectory (excluding artificial links).}
#' \item{`Max_density`}{Maximum number of segments represented by at least one of
#' the segments of the representative trajectory (excluding artificial links).}
#' \item{`Avg_depth`}{Mean value of the k-d tree depths, that is, the number of
#' partitions of the ordination space until finding a region with `minSegs` segments
#' or less.}
#' \item{`Max_depth`}{Maximum depth in the k-d tree, that is, the number of partitions
#' of the ordination space until finding a region with `minSegs` segments or less.}
#' }
#'
#' @seealso
#' [retra_edr()] for identifying representative trajectories in EDRs applying
#' RETRA-EDR.
#'
#' [`define_retra()`] for generating an object of class `RETRA` from trajectory
#' features.
#'
#'
#' @export
#'
#' @examples
#' # Apply RETRA-EDR to identify representative trajectories
#' d = EDR_data$EDR1$state_dissim
#' trajectories = EDR_data$EDR1$abundance$traj
#' states = EDR_data$EDR1$abundance$state
#' RT <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)
#'
#' # Summarize the properties of the representative trajectories in a data frame
#' summary(RT)
#'
summary.RETRA <- function(object, ...) {
  df <- data.frame(ID = names(object),
                   Size = vapply(object, function(x){x$Size}, numeric(1)),
                   Length = vapply(object, function(x){x$Length}, numeric(1)),
                   Avg_link = vapply(object, function(x) {
                     if (all(!is.na(x$Link_distance))) {
                       mean(x$Link_distance$Distance)
                     } else {NA}
                   }, numeric(1)),
                   Sum_link = vapply(object, function(x) {
                     if (all(!is.na(x$Link_distance))) {
                       sum(x$Link_distance$Distance)
                     } else {NA}
                   }, numeric(1)),
                   Avg_density = vapply(object, function(x) {
                     if (all(!is.na(x$Seg_density))) {
                       mean(x$Seg_density$Density)
                     } else {NA}
                   }, numeric(1)),
                   Max_density = vapply(object, function(x) {
                     if (all(!is.na(x$Seg_density))) {
                       max(x$Seg_density$Density)
                     } else {NA}
                   }, numeric(1)),
                   Avg_depth = vapply(object, function(x) {
                     if (all(!is.na(x$Seg_density))) {
                       mean(x$Seg_density$kdTree_depth)
                     } else {NA}
                   }, numeric(1)),
                   Max_depth = vapply(object, function(x) {
                     if (all(!is.na(x$Seg_density))) {
                       max(x$Seg_density$kdTree_depth)
                     } else {NA}
                   }, numeric(1)))
  return(df)
}
