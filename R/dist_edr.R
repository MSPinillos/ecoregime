#' Dissimilarities between Ecological Dynamic Regimes
#'
#' @description
#' Generate a matrix containing dissimilarities between one or more pairs of
#' Ecological Dynamic Regimes (EDR). `dist_edr()` computes different dissimilarity
#' indices, all of them based on the dissimilarities between the trajectories of
#' two EDRs.
#'
#' @param d Symmetric matrix or object of class [`dist`] containing the
#' dissimilarities between each pair of states of all trajectories in the EDR or
#' the dissimilarities between each pair of trajectories.
#' @param d.type One of `"dStates"` (if `d` contains state dissimilarities) or
#' `"dTraj"` (if `d` contains trajectory dissimilarities).
#' @param trajectories Only if `d.type` = `"dStates"`. Vector indicating the
#' trajectory or site corresponding to each entry in `d`.
#' @param states Only if `d.type` = `"dStates"`. Vector of integers indicating the
#' order of the states in `d` for each trajectory.
#' @param edr Vector indicating the EDR to which each trajectory/state in `d` belongs.
#' @param metric A string indicating the dissimilarity index to be used: `"dDR"`
#' (default), `"minDist"`, `"maxDist"`.
#' @param symmetrize String naming the function to be called to symmetrize the
#' resulting dissimilarity matrix (`"mean"`, `"min"`, `"max`, `"lower"`, `"upper"`).
#' If `NULL` (default), the matrix is not symmetrized.
#' @param ... Only if `d.type` = `"dStates"`. Further arguments to calculate
#' trajectory dissimilarities. See [ecotraj::trajectoryDistances()].
#'
#' @details
#' The implemented metrics are:
#'
#' \describe{
#' \item{`"dDR"`}{
#' \eqn{
#' d_{DR}(R_1, R_2) = \frac{1}{n} \sum_{i=1}^{n} d_{TR}(T_{1i}, R_2)
#' }
#' }
#'
#' \item{`"minDist"`}{
#' \eqn{
#' d_{DRmin}(R_1, R_2) = \min_{i=1}^{n} \{ d_{TR}(T_{1i}, R_2) \}
#' }
#' }
#'
#' \item{`"maxDist"`}{
#' \eqn{
#' d_{DRmax}(R_1, R_2) = \max_{i=1}^{n} \{ d_{TR}(T_{1i}, R_2) \}
#' }
#' }
#' }
#'
#' where \eqn{R_1} and \eqn{R_2} are two EDRs composed of \eqn{n} and \eqn{m}
#' ecological trajectories, respectively, and \eqn{d_{TR}(T_{1i}, R_2)} is the
#' dissimilarity between the trajectory \eqn{T_{1i}} of \eqn{R_1} and the closest
#' trajectory of \eqn{R_2}:
#'
#' \eqn{
#' d_{TR}(T_{1i}, R_2) = \min\{d_T(T_{1i}, T_{21}), ... , d_T(T_{1i}, T_{2m})\}
#' }
#'
#' The metrics calculated are not necessarily symmetric. That is, \eqn{d_{DR}(R_1, R_2)}
#' is not necessarily equal to \eqn{d_{DR}(R_2, R_1)}. It is possible to symmetrize
#' the returned matrix by indicating the name of the function to be used in `symmetrize`:
#'
#' \describe{
#' \item{`"mean"`}{
#' \eqn{
#' d_{DRsym} = \frac{d_{DR}(R_1, R_2) + d_{DR}(R_2, R_1)}{2}
#' }
#' }
#'
#' \item{`"min"`}{
#' \eqn{
#' d_{DRsym} = \min\{d_{DR}(R_1, R_2), d_{DR}(R_2, R_1)\}
#' }
#' }
#'
#' \item{`"max"`}{
#' \eqn{
#' d_{DRsym} = \max\{d_{DR}(R_1, R_2), d_{DR}(R_2, R_1)\}
#' }
#' }
#'
#' \item{`"lower"`}{
#' The lower triangular part of the dissimilarity matrix is used.
#' }
#'
#' \item{`"upper"`}{
#' The upper triangular part of the dissimilarity matrix is used.
#' }
#'
#' }
#'
#' @return
#' Matrix including the dissimilarities between every pair of EDRs.
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023. Ecological Dynamic
#' Regimes: Identification, characterization, and comparison. *Ecological Monographs*.
#' <doi:10.1002/ecm.1589>
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
#' dStates <- vegan::vegdist(abun[, -c(1, 2, 3)], method = "bray")
#'
#' # Use the labels in dStates to define the trajectories to which each state
#' # belongs
#' id_traj <- vapply(strsplit(labels(dStates), "_"), function(x){
#'                     paste0(x[1], "_", x[2])
#'                 }, character(1))
#' id_state <- vapply(strsplit(labels(dStates), "_"), function(x){
#'                     as.integer(x[3])
#'                 }, integer(1))
#' id_edr <- vapply(strsplit(labels(dStates), "_"), function(x){
#'                     paste0("EDR", x[1])
#'                 }, character(1))
#'
#' # Calculate dissimilarities between every pair of trajectories
#' dTraj <- ecotraj::trajectoryDistances(d = dStates, sites = id_traj,
#'                                       surveys = id_state,
#'                                       distance.type = "DSPD")
#'
#' # Use labels in dTraj to identify EDRs
#' id_edr_traj <- vapply(strsplit(labels(dTraj), "_"), function(x){
#'                     paste0("EDR", x[1])
#'                 }, character(1))
#'
#' # Compute dissimilarities between EDRs:
#' # 1.1) without symmetrizing the matrix and using state dissimilarities
#' dEDR <- dist_edr(d = dStates, d.type = "dStates",
#'                  trajectories = id_traj, states = id_state, edr = id_edr,
#'                  metric = "dDR", symmetrize = NULL)
#'
#' # 1.2) without symmetrizing the matrix and using trajectory dissimilarities
#' dEDR <- dist_edr(d = dTraj, d.type = "dTraj", edr = id_edr_traj,
#'                  metric = "dDR", symmetrize = NULL)
#'
#' # 2) symmetrizing by averaging elements on and below the diagonal
#' dEDR <- dist_edr(d = dTraj, d.type = "dTraj", edr = id_edr_traj,
#'                  metric = "dDR", symmetrize = "mean")
#' # 3) symmetrizing by using the minimum dissimilarity between two EDRs
#' dEDR <- dist_edr(d = dTraj, d.type = "dTraj", edr = id_edr_traj,
#'                  metric = "dDR", symmetrize = "min")
#'
dist_edr <- function(d, d.type, trajectories = NULL, states = NULL, edr, metric = "dDR", symmetrize = NULL, ...){

  metric <- match.arg(metric, c("dDR", "minDist", "maxDist"))
  if (!is.null(symmetrize)) {
    symmetrize <- match.arg(symmetrize, c("mean", "min", "max", "lower", "upper"))
  }

  d.type <- match.arg(d.type, c("dStates", "dTraj"))
  if (d.type == "dTraj") {
    states <- NULL
  }

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for d, trajectories, states, and reference
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (d.type == "dStates") {
    if (is.null(trajectories)) {
      stop("If 'd.type' = \"dStates\", you must provide a value for 'trajectories'.")
    }
    if (is.null(states)) {
      stop("If 'd.type' = \"dStates\", you must provide a value for 'states'.")
    }
    if (length(trajectories) != nrow(as.matrix(d))) {
      stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
    }
    if (length(states) != nrow(as.matrix(d))) {
      stop("The length of 'states' must be equal to both dimensions in 'd'.")
    }
  }

  # Check input data is correct
  if(length(edr) != nrow(as.matrix(d))){
    stop("'edr' needs to have a length equal to the number of rows and columns of 'd'. \nProvide the EDR of all trajectories considered in 'd'.")
  }

  ## TRAJECTORIES & DISSIMILARITIES --------------------------------------------

  # Trajectory dissimilarity
  if(d.type == "dStates"){
    edr.df <- unique(data.frame(traj = paste0(edr, "_", trajectories), edr = edr))
    edr <- edr.df$edr
    dTrajmat <- as.matrix(ecotraj::trajectoryDistances(d = d,
                                                       sites = paste0(edr, "_", trajectories),
                                                       surveys = states,...))
  }
  if(d.type == "dTraj") {
    dTrajmat <- as.matrix(d)
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
      if(metric == "dDR"){
        dReg[iR1, iR2] <- sum(minD_T1R2)/(Ntraj_edr[iR1])           # dDR(R1, R2)
      }
      if(metric == "minDist"){
        dReg[iR1, iR2] <- min(minD_T1R2)
      }
      if(metric == "maxDist"){
        dReg[iR1, iR2] <- max(minD_T1R2)
      }
    }
  }

  # Symmetrize dDR if required
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
