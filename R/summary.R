#' Summarize representative trajectories
#'
#' @description
#' Summary of the properties of each representative trajectory obtained
#' with RETRA-EDR.
#'
#' @param object Object of class `RETRA` returned from function [retra_edr()].
#' @param ... (not used)
#'
#' @return
#' Data frame with nine columns and one row for each representative trajectory
#' in `object`. The columns in the returned data frame contain the following
#' information:
#' * `ID`: Identifier of the representative trajectories returned from [retra_edr()].
#' * `Size`: Number of states forming each representative trajectory.
#' * `Length`: Sum of the dissimilarities in `d` between every pair of consecutive
#' states forming the representative trajectories.
#' * `Avg_link`: Mean value of the dissimilarities between every pair of states
#' linking two representative segments in the representative trajectories.
#' * `Sum_link`: Sum of the dissimilarities between every pair of states linking
#' two representative segments in the representative trajectories.
#' * `Avg_density`: Mean value of the number of segments represented by each real
#' segment of the representative trajectory (i.e., excluding links).
#' * `Max_density`: Maximum number of segments represented by at least of the
#' segments of the representative trajectory (excluding links).
#' * `Avg_depth`: Mean value of the kd-tree depths, that is, the number of partitions
#' of the ordination space until finding a region with `minSegs` segments or less.
#' * `Max_depth`: Maximum depth in the kd-tree, that is, the number of partitions
#' of the ordination space until finding a region with `minSegs` segments or less.
#'
#' @note
#' Note that the `density` of the dense regions identified in the EDR tend to `minSegs`.
#' In contrast, greater values of `kdTree_depth` indicate a higher stability of
#' dense regions over sequential partitions of the ordination space.
#'
#' @seealso
#' [retra_edr()] for identifying representative trajectories in EDRs.
#'
#' @export
#'
#' @examples
#' d = EDR_data$EDR1$state_dissim
#' trajectories = EDR_data$EDR1$abundance$traj
#' states = EDR_data$EDR1$abundance$state
#' RT <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)
#' summary(RT)
#'
summary.RETRA <- function(object, ...) {
  df <- data.frame(ID = names(object),
                   Size = vapply(object, function(x){x$Size}, numeric(1)),
                   Length = vapply(object, function(x){x$Length}, numeric(1)),
                   Avg_link = vapply(object, function(x) {
                     mean(x$Link_distance$Distance)
                   }, numeric(1)),
                   Sum_link = vapply(object, function(x) {
                     sum(x$Link_distance$Distance)
                   }, numeric(1)),
                   Avg_density = vapply(object, function(x) {
                     mean(x$Seg_density$Density)
                   }, numeric(1)),
                   Max_density = vapply(object, function(x) {
                     max(x$Seg_density$Density)
                   }, numeric(1)),
                   Avg_depth = vapply(object, function(x) {
                     mean(x$Seg_density$kdTree_depth)
                   }, numeric(1)),
                   Max_depth = vapply(object, function(x) {
                     max(x$Seg_density$kdTree_depth)
                   }, numeric(1)))
  return(df)
}
