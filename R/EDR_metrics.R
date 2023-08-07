#' Metrics of Ecological Dynamic Regimes
#'
#' @description
#' Set of metrics to analyze the distribution and variability of trajectories
#' in Ecological Dynamic Regimes (EDR), including dynamic dispersion (dDis), dynamic
#' beta diversity (dBD), and dynamic evenness (dEve).
#'
#' @name EDR_metrics
#' @aliases dDis dBD dEve
#'
#' @param d Symmetric matrix or object of class [`dist`] containing the
#' dissimilarities between each pair of states of all trajectories in the EDR or
#' the dissimilarities between each pair of trajectories. To compute dDis, `d`
#' needs to include the dissimilarities between all states/trajectories and the
#' states/trajectory of reference.
#' @param d.type One of `"dStates"` (if `d` contains state dissimilarities) or
#' `"dTraj"` (if `d` contains trajectory dissimilarities).
#' @param trajectories Vector indicating the trajectory or site corresponding to
#' each entry in `d`.
#' @param states Only if `d.type` = `"dStates"`. Vector of integers indicating the
#' order of the states in `d` for each trajectory.
#' @param reference Vector of the same class as `trajectories` and length equal
#' to one, indicating the reference trajectory to compute dDis.
#' @param w.type Method used to weight individual trajectories:
#' * `"none"`: All trajectories are considered equally relevant (default).
#' * `"length"`: Trajectories are weighted by their length, calculated as the
#' sum of the dissimilarities between every pair of consecutive states. `d` must
#' contain dissimilarities between trajectory states and `d.type` = `"dStates"`.
#' * `"size"`: Trajectories are weighted by their size, calculated as the number
#' of states forming the trajectory. `d` must contain dissimilarities between
#' trajectory states and `d.type` = `"dStates"`.
#' * `"precomputed"`: Trajectories weighted according to different criteria.
#' @param w.values Only if `w.type` = `"precomputed"`. Numeric vector of length
#' equal to the number of trajectories containing the weight of each trajectory.
#' @param ... Only if `d.type` = `"dStates"`. Further arguments to calculate
#' trajectory dissimilarities. See [ecotraj::trajectoryDistances()].
#'
#' @return
#' * `dDis()` returns the value of dynamic dispersion for a given trajectory
#' taken as a reference.
#' * `dBD()` returns the value of dynamic beta diversity.
#' * `dEve()` returns the value of dynamic evenness.
#'
#' @details
#'
#' \strong{Dynamic dispersion (`dDis()`)}
#'
#' *dDis* is calculated as the average dissimilarity between each trajectory in an
#' EDR and a target trajectory taken as reference (Sánchez-Pinillos et al., 2023).
#'
#' \eqn{
#' dDis = \frac{\sum_{i=1}^{m}d_{i\alpha}}{m}
#' }
#'
#' where \eqn{d_{i\alpha}} is the dissimilarity between trajectory \eqn{i} and
#' the trajectory of reference \eqn{\alpha}, and \eqn{m} is the number of trajectories.
#'
#' Alternatively, it is possible to calculate a weighted mean of the dissimilarities
#' by assigning a weight to each trajectory.
#'
#' \eqn{
#' dDis = \frac{\sum_{i=1}^{m}w_{i}d_{i\alpha}}{\sum_{i=1}^{m}w_{i}}
#' }
#'
#' where \eqn{w_{i}} is the weight assigned to trajectory \eqn{i}.
#'
#'
#' \strong{Dynamic beta diversity (`dBD()`)}
#'
#' *dBD* quantifies the overall variation of the trajectories in an EDR and is
#' equivalent to the average distance to the centroid of the EDR
#' (De Cáceres et al., 2019).
#'
#' \eqn{
#' dBD = \frac{\sum_{i=1}^{m-1}\sum_{j=i+1}^{m}d_{ij}^{2}}{m(m-1)}
#' }
#'
#' \strong{Dynamic evenness (`dEve()`)}
#'
#' *dEve* quantifies the regularity with which an EDR is filled by the individual
#' trajectories (Sánchez-Pinillos et al., 2023).
#'
#' \eqn{
#' dEve = \frac{\sum_{l=1}^{m-1}\min(\frac{d_{ij}}{\sum_{l=1}^{m-1}d_{ij}}, \frac{1}{m-1}) - \frac{1}{m-1}}{1-\frac{1}{1-1}}
#' }
#'
#' where \eqn{d_{ij}} is the dissimilarity between trajectories \eqn{i} and \eqn{j}
#' linked in a minimum spanning tree by the link \eqn{l}.
#'
#' Optionally, it is possible to weight the trajectories of the EDR. In that case,
#' *dEve* becomes analogous to the functional evenness index proposed by Villéger
#' et al. (2008).
#'
#' \eqn{
#' dEve_{w} = \frac{\sum_{l=1}^{m-1}\min(\frac{EW_{ij}}{\sum_{l=1}^{m-1}EW_{ij}}, \frac{1}{m-1}) - \frac{1}{m-1}}{1-\frac{1}{1-1}}
#' }
#'
#' where \eqn{EW_{ij}} is the weighted evenness:
#'
#' \eqn{
#' EW_{ij} = \frac{d_{ij}}{w_i + w_j}
#' }
#'
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' De Cáceres, M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit R &
#' Hubbell S. (2019). Trajectory analysis in community ecology. Ecological
#' Monographs.
#'
#' Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023. Ecological Dynamic
#' Regimes: Identification, characterization, and comparison. *Ecological Monographs*.
#' <doi:10.1002/ecm.1589>
#'
#' Villéger, S., Mason, N.W.H., Mouillot, D. (2008) New multidimensional functional
#' diversity indices for a multifaced framework in functional ecology. Ecology.
#'
#' @export
#'
#' @examples
#' # Data to compute dDis, dBD, and dEve
#' dStates <- EDR_data$EDR1$state_dissim
#' dTraj <- EDR_data$EDR1$traj_dissim
#' trajectories <- paste0("T", EDR_data$EDR1$abundance$traj)
#' states <- EDR_data$EDR1$abundance$state
#'
#' # Dynamic dispersion taking the first trajectory as reference
#' dDis(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories),
#'          reference = "T1")
#'
#' # Dynamic dispersion weighting trajectories by their length
#' dDis(d = dStates, d.type = "dStates", trajectories = trajectories, states = states,
#'          reference = "T1", w.type = "length")
#'
#' # Dynamic beta diversity using trajectory dissimilarities
#' dBD(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories))
#'
#' # Dynamic evenness
#' dEve(d = dStates, d.type = "dStates", trajectories = trajectories, states = states)
#'
#' # Dynamic evenness considering that the 10 first trajectories are three times
#' # more relevant than the rest
#' w.values <- c(rep(3, 10), rep(1, length(unique(trajectories))-10))
#' dEve(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories),
#'          w.type = "precomputed", w.values = w.values)
#'


#### DYNAMIC DISPERSION (dDis) ####

dDis <- function(d, d.type, trajectories, states = NULL, reference, w.type = "none", w.values, ...){

  d.type <- match.arg(d.type, c("dStates", "dTraj"))
  w.type <- match.arg(w.type, c("none", "length", "size", "precomputed"))
  if (d.type == "dTraj") {
    states <- NULL
  }

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for d, trajectories, states, and reference
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }

  if (d.type == "dStates") {
    if (is.null(states)) {
      stop("If 'd.type' = \"dStates\", you must provide a value for 'states'.")
    }
    if (length(states) != nrow(as.matrix(d))) {
      stop("The length of 'states' must be equal to both dimensions in 'd'.")
    }
  }

  if (length(reference) != 1) {
    stop("'reference' needs to have a length equal to one.")
  }
  if (sum(reference %in% trajectories) == 0) {
    if (d.type == "dStates") {
      stop("'reference' needs to be specified in 'trajectories' and 'd' must include dissimilarities between the observations of the reference trajectory and those of the other trajectories.")
    }
    if (d.type == "dTraj") {
      stop("'reference' needs to be specified in 'trajectories' and 'd' must include dissimilarities between the reference trajectory and the other trajectories.")
    }
  }

  ## INDIVIDUAL TRAJECTORIES ---------------------------------------------------

  # Identify individual trajectories (no ref)
  noRef <- trajectories[trajectories != reference]
  iRef = which(unique(trajectories) == reference)
  if(!is.null(states)){
    noRef_states <- states[trajectories != reference]
  }

  unique_noRef <- unique(noRef)
  Ntraj <- length(unique_noRef)

  ## WEIGHTING TRAJECTORIES ----------------------------------------------------

  # Check the format for w.type and w.values
  if (w.type == "precomputed") {
    if (length(w.values) == 1) {
      warning("'w.values' has length 1. Equal weights will be assigned to all trajectories.")
      w.type = "none"
    } else if (length(w.values) != Ntraj){
      stop("The length of 'w.values' needs to be equal to the number of trajectories to be evaluated (excluding the reference trajectory).")
    }
    if(!is.numeric(w.values)){
      stop("'w.values' needs to be numeric")
    }
  }

  if (w.type == "none") {
    w.values <- rep(1, Ntraj)
  }

  if (w.type == "length") {
    if (d.type == "dTraj") {
      stop("If w.type = \"length\", 'd' needs to contain dissimilarities between trajectory states (i.e., d.type =\"dStates\").")
    } else {
      trajL = ecotraj::trajectoryLengths(d = as.matrix(d)[trajectories %in% noRef, trajectories %in% noRef],
                                         sites = noRef, surveys = noRef_states)
      w.values = trajL$Trajectory
    }
  }

  if (w.type == "size") {
    if (d.type == "dTraj") {
      stop("If w.type = \"size\", 'd' needs to contain dissimilarities between trajectory states (i.e., d.type =\"dStates\").")
    } else {
      w.values = table(noRef)
    }
  }

  ## COMPUTE dDis --------------------------------------------------------------

  # If d = "dStates", calculate trajectory dissimilarities
  if (d.type == "dStates") {
    trajD.Ref = as.matrix(ecotraj::trajectoryDistances(d = d, sites = trajectories,
                                                       surveys = states,...))
  } else if (d.type == "dTraj") {
    trajD.Ref = as.matrix(d)
    dimnames(trajD.Ref) <- list(trajectories, trajectories)
  }

  # Compute dDis
  dtoRef = trajD.Ref[reference, ][-iRef]
  dDis = weighted.mean(dtoRef, w = w.values)
  names(dDis) <- paste0("dDis (ref. ", reference, ")")

  return(dDis)

}

################################################################################

#' @rdname EDR_metrics
#' @export

#### DYNAMIC BETA DIVERSITY (dBD) ####

dBD <- function(d, d.type, trajectories, states = NULL, ...){

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

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }

  if (d.type == "dStates") {
    if (is.null(states)) {
      stop("If 'd.type' = \"dStates\", you must provide a value for 'states'.")
    }
    if (length(states) != nrow(as.matrix(d))) {
      stop("The length of 'states' must be equal to both dimensions in 'd'.")
    }
  }

  ## TRAJECTORIES & DISSIMILARITIES --------------------------------------------

  # Number of trajectories
  unique_traj <- unique(trajectories)
  Ntraj <- length(unique_traj)

  # Trajectory dissimilarity
  if(d.type == "dStates"){
    trajD <- as.dist(ecotraj::trajectoryDistances(d = d, sites = trajectories, surveys = states,...))
  }
  if(d.type == "dTraj") {
    trajD <- as.dist(d)
  }

  ## COMPUTE dBD ---------------------------------------------------------------

  # Compute dBD
  SS <- (sum(trajD^2)) / Ntraj
  dBD <- SS / (Ntraj - 1)
  names(dBD) <- "dBD"

  return(dBD)

}

################################################################################

#' @rdname EDR_metrics
#' @export

#### DYNAMIC EVENNESS (dEve) ####

dEve <- function(d, d.type, trajectories, states = NULL, w.type = "none", w.values, ...){

  d.type <- match.arg(d.type, c("dStates", "dTraj"))
  w.type <- match.arg(w.type, c("none", "length", "size", "precomputed"))
  if (d.type == "dTraj") {
    states <- NULL
  }

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for d, trajectories, states, and reference
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }

  if (d.type == "dStates") {
    if (is.null(states)) {
      stop("If 'd.type' = \"dStates\", you must provide a value for 'states'.")
    }
    if (length(states) != nrow(as.matrix(d))) {
      stop("The length of 'states' must be equal to both dimensions in 'd'.")
    }
  }

  ## TRAJECTORIES & DISSIMILARITIES --------------------------------------------

  # Number of trajectories
  unique_traj <- unique(trajectories)
  Ntraj <- length(unique_traj)

  # Trajectory dissimilarity
  if(d.type == "dStates"){
    trajD <- as.matrix(ecotraj::trajectoryDistances(d = as.matrix(d), sites = trajectories, surveys = states,...))
  }
  if(d.type == "dTraj") {
    trajD <- as.matrix(d)
  }

  ## WEIGHTING TRAJECTORIES ----------------------------------------------------

  # Check the format for w.type and w.values
  if (w.type == "precomputed") {
    if (length(w.values) == 1) {
      warning("'w.values' has length 1. Equal weights will be assigned to all trajectories.")
      w.type = "none"
    } else if (length(w.values) != Ntraj){
      stop("The length of 'w.values' needs to be equal to the number of trajectories to be evaluated (excluding the reference trajectory).")
    }
    if(!is.numeric(w.values)){
      stop("'w.values' needs to be numeric")
    }
  }

  if (w.type == "none") {
    w.values <- rep(1, Ntraj)
  }

  if (w.type == "length") {
    if (d.type == "dTraj") {
      stop("If w.type = \"length\", 'd' needs to contain dissimilarities between trajectory states (i.e., d.type =\"dStates\").")
    } else {
      trajL = ecotraj::trajectoryLengths(d = d, sites = trajectories, surveys = states)
      w.values = trajL$Trajectory
    }
  }

  if (w.type == "size") {
    if (d.type == "dTraj") {
      stop("If w.type = \"size\", 'd' needs to contain dissimilarities between trajectory states (i.e., d.type =\"dStates\").")
    } else {
      w.values = table(trajectories)
    }
  }

  ## COMPUTE dEve ---------------------------------------------------------------

  # Identify trajectories connected in a MST
  i <-  2
  n <- Ntraj
  row = numeric()
  col = numeric()
  while (i <= Ntraj) {
    row = c(row, i:Ntraj)
    col = c(col, rep(i-1, n-1))
    i = i + 1
    n = n - 1
  }
  MST.links <- data.frame(row = row, col = col, linked = as.numeric(as.dist(ape::mst(trajD))))
  positive.links <- MST.links[which(MST.links$linked > 0), ]

  # Calculate the weighted evenness (EW) and the partial weighted evenness (PEW)
  dlinks <- numeric()
  Sum_wi_wj <- numeric()
  for(ilink in 1:(Ntraj-1)){
    dlinks <- c(dlinks, trajD[positive.links$row[ilink], positive.links$col[ilink]])
    Sum_wi_wj <- c(Sum_wi_wj, sum(w.values[c(positive.links$row[ilink], positive.links$col[ilink])]))
  }
  EW <- dlinks / Sum_wi_wj
  PEW <- EW / sum(EW)

  # Compute dEve
  min_PEW <- vapply(PEW, function(iPEW){
    min(iPEW, 1 / (Ntraj - 1))
  }, numeric(1))
  dEve <- (sum(min_PEW) - 1/(Ntraj-1)) / (1 - 1/(Ntraj-1))
  names(dEve) <- "dEve"

  return(dEve)

}





