#' Dissimilarities between Ecological Dynamic Regimes
#'
#' @description
#' Generate a matrix compiling dissimilarities between one or more pairs of
#' Ecological Dynamic Regimes (EDR). `dist_edr` computes different dissimilarity
#' indices, all of them based on the dissimilarities between the trajectories of
#' two EDRs.
#'
#' @param dTraj Either a square matrix or an object of class [`dist`] containing
#' the dissimilarities between every pair of trajectories in the target EDRs.
#' @param edr Vector indicating the EDR to which each trajectory in `dTraj` belongs.
#' @param metric A string indicating the dissimilarity index to be used: `"DDR"`
#' (default), `"minDist"`, `"maxDist"`.
#' @param symmetrize String naming the function to be called to symmetrize the
#' resulting dissimilarity matrix (`"mean"`, `"min"`, `"max`, `"lower"`, `"upper"`).
#' If `NULL` (default), the matrix is not symmetrized.
#'
#' @details
#' The implemented metrics are:
#' * `"DDR"`: \eqn{
#' D_{DR}(R_1, R_2) = \frac{1}{n} \sum_{i=1}^{n} D_{TR}(T_{1i}, R_2)
#' }
#' * `"minDist"`: \eqn{
#' D_{DRmin}(R_1, R_2) = \min_{i=1}^{n} \{ D_{TR}(T_{1i}, R_2) \}
#' }
#' * `"maxDist"`: \eqn{
#' D_{DRmax}(R_1, R_2) = \max_{i=1}^{n} \{ D_{TR}(T_{1i}, R_2) \}
#' }
#'
#' where \eqn{R_1} and \eqn{R_2} are two EDRs composed of \eqn{n} and \eqn{m}
#' ecological trajectories, respectively, and \eqn{D_{TR}(T_{1i}, R_2)} is the
#' dissimilarity between the trajectory \eqn{T_{1i}} of \eqn{R_1} and the closest
#' trajectory of \eqn{R_2} (Sánchez-Pinillos et al.):
#'
#' \eqn{
#' D_{TR}(T_{1i}, R_2) = \min\{D_T(T_{1i}, T_{21}), ... , D_T(T_{1i}, T_{2m})\}
#' }
#'
#' The metrics calculated are not necessarily symmetric. That is, \eqn{D_{DR}(R_1, R_2)}
#' is not necessarily equal to \eqn{D_{DR}(R_2, R_1)}. It is possible to symmetrize
#' the returned matrix by indicating the name of the function to be used in `symmetrize`:
#'
#' * `"mean"`: \eqn{
#' D_{DRsym} = \frac{D_{DR}(R_1, R_2) + D_{DR}(R_2, R_1)}{2}
#' }
#' * `"min"`: \eqn{
#' D_{DRsym} = \min\{D_{DR}(R_1, R_2), D_{DR}(R_2, R_1)\}
#' }
#' * `"max"`: \eqn{
#' D_{DRsym} = \max\{D_{DR}(R_1, R_2), D_{DR}(R_2, R_1)\}
#' }
#' * `"lower"`: The lower triangular part of the dissimilarity matrix is used.
#' * `"upper"`: The upper triangular part of the dissimilarity matrix is used.
#'
#'
#' @return
#' Matrix including the dissimilarities between every pair of EDRs.
#'
#' @author Martina Sánchez-Pinillos, CNRS, Univ. Montpellier
#'
#' @references
#' Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. Ecological Dynamic
#' Regimes: Identification, characterization, and comparison
#'
#'
#' @export
#'
#' @examples
#' # Load species abundances and compile in a data frame
#' abun1 <- EDR_data$EDR1$abundance
#' abun2 <- EDR_data$EDR2$abundance
#' abun3 <- EDR_data$EDR3$abundance
#' abun <- data.frame(rbind(abun1, abun2, abun3))
#'
#' # Define row names in abun to keep the reference of the EDR, trajectory, and
#' # state
#' row.names(abun) <- paste0(abun$EDR, "_", abun$traj, "_", abun$state)
#'
#' # Calculate dissimilarities between every pair of states
#' # For example, Bray-Curtis index
#' library(vegan)
#' dStates <- vegdist(abun[, -c(1, 2, 3)], method = "bray")
#'
#' # Use the labels in dStates to define the trajectories to which each state
#' # belongs
#' id_traj <- vapply(strsplit(labels(dStates), "_"), function(x){
#'                     paste0(x[1], "_", x[2])
#'                 }, character(1))
#' id_state <- vapply(strsplit(labels(dStates), "_"), function(x){
#'                     as.integer(x[3])
#'                 }, integer(1))
#'
#' # Calculate dissimilarities between every pair of trajectories
#' library(ecotraj)
#' dTraj <- trajectoryDistances(d = dStates,
#'                              sites = id_traj,
#'                              surveys = id_state,
#'                              distance.type = "DSPD")
#'
#' # Use labels in dTraj to identify EDRs
#' id_edr <- vapply(strsplit(labels(dTraj), "_"), function(x){
#'                     paste0("EDR", x[1])
#'                 }, character(1))
#'
#' # Compute dissimilarities between EDRs
#' # without symmetrizing the matrix
#' dEDR <- dist_edr(dTraj = dTraj,
#'                  edr = id_edr,
#'                  metric = "DDR",
#'                  symmetrize = NULL)
#' # symmetrizing by averaging elements on and below the diagonal
#' dEDR <- dist_edr(dTraj = dTraj,
#'                  edr = id_edr,
#'                  metric = "DDR",
#'                  symmetrize = "mean")
#' # symmetrizing by using the minimum dissimilarity between two EDRs
#' dEDR <- dist_edr(dTraj = dTraj,
#'                  edr = id_edr,
#'                  metric = "DDR",
#'                  symmetrize = "min")
#' # symmetrizing by using the lower triangular part of the asymmetric matrix
#' dEDR <- dist_edr(dTraj = dTraj,
#'                  edr = id_edr,
#'                  metric = "DDR",
#'                  symmetrize = "lower")
#'
dist_edr <- function(dTraj, edr, metric = "DDR", symmetrize = NULL){

  metric <- match.arg(metric, c("DDR", "minDist", "maxDist"))
  if (!is.null(symmetrize)) {
    symmetrize <- match.arg(symmetrize, c("mean", "min", "max", "lower", "upper"))
  }

  # Convert dTraj into a matrix
  if (all(!is.matrix(dTraj), !dendextend::is.dist(dTraj)) |
      dim(dTraj)[1] != dim(dTraj)[2]) {
    stop("'d' must be a square matrix or an object of class 'dist'.")
  }
  dTrajmat <- as.matrix(dTraj)

  # Check input data is correct
  if(length(edr) != nrow(dTrajmat)){
    stop(cat("'edr' needs to have a length equal to the number of rows and columns of 'dTraj'. \nProvide the EDR of each trajectory considered in 'dTraj'.", "\n"))
  }

  # Nb edr and trajectories/edr
  ID_edr <- unique(edr)
  Ntraj_edr <- table(edr)

  # Empty matrix to compile regime distances
  dReg <- matrix(0, ncol = length(ID_edr), nrow = length(ID_edr), dimnames = list(ID_edr, ID_edr))


  for (iR1 in seq_along(ID_edr)) {                      # R1

    for (iR2 in seq_along(ID_edr)) {                    # R2
      minD_T1R2 <- numeric()

      for (T1 in 1:Ntraj_edr[iR1]) {        # for each T1i, find the minimum distance to R2; min{DSDSP(T1i, T21), ..., DSDSP(T1i, T2m)}
        iT1 <- sum(Ntraj_edr[1:iR1]) - Ntraj_edr[iR1] + T1                         # index of T1i in dTraj
        T2 <- (sum(Ntraj_edr[1:iR2]) - Ntraj_edr[iR2] + 1):sum(Ntraj_edr[1:iR2])  # indices of all T21,..., T2m in dTraj
        minD_T1R2 <- c(minD_T1R2, min(dTrajmat[iT1, T2]))            # min{DSDSP(T1i, T21), ..., DSDSP(T1i, T2m)}
      }
      if(metric == "DDR"){
        dReg[iR1, iR2] <- sum(minD_T1R2)/(Ntraj_edr[iR1])           # DDR(R1, R2)
      }
      if(metric == "minDist"){
        dReg[iR1, iR2] <- min(minD_T1R2)
      }
      if(metric == "maxDist"){
        dReg[iR1, iR2] <- max(minD_T1R2)
      }
    }
  }

  # Symmetrize DDR if required
  if (!is.null(symmetrize)) {
    if (symmetrize == "mean") {
      for (iR1 in seq_along(ID_edr)) {
        for (iR2 in seq_along(ID_edr)) {
          dReg[iR1, iR2] <- mean(c(dReg[iR1, iR2], dReg[iR2, iR1]))
          dReg[iR2, iR1] <- dReg[iR1, iR2]
        }
      }
    }
    if (symmetrize == "min") {
      for (iR1 in seq_along(ID_edr)) {
        for (iR2 in seq_along(ID_edr)) {
          dReg[iR1, iR2] <- min(c(dReg[iR1, iR2], dReg[iR2, iR1]))
          dReg[iR2, iR1] <- dReg[iR1, iR2]
        }
      }
    }
    if (symmetrize == "max") {
      for (iR1 in seq_along(ID_edr)) {
        for (iR2 in seq_along(ID_edr)) {
          dReg[iR1, iR2] <- max(c(dReg[iR1, iR2], dReg[iR2, iR1]))
          dReg[iR2, iR1] <- dReg[iR1, iR2]
        }
      }
    }
    if (symmetrize == "lower") {
      dReg[upper.tri(dReg)] = t(dReg)[upper.tri(dReg)]
    }
    if (symmetrize == "upper") {
      dReg[lower.tri(dReg)] = t(dReg)[lower.tri(dReg)]
    }
  }

  return(dReg)

}
