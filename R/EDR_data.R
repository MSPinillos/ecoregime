#' Ecological Dynamic Regime data
#'
#' @description
#' Example datasets to characterize and compare EDRs, including abundance data;
#' state, segment, and trajectory dissimilarity matrices for 90 artificial communities
#' belonging to three different EDRs.
#'
#' @format
#' List of three nested sublists (`"EDR1"`, `"EDR2"`, and `"EDR3"`), each
#' associated with one EDR, including the following elements:
#' * `abundance`: Data table with 15 columns and one row for each community state:
#'     + `EDR`: Integer indicating the identifier of the EDR.
#'     + `traj`: Integer containing the identifier of the trajectory for each
#'     artificial community in the corresponding EDR. Each trajectory represents
#'     a different sampling unit.
#'     + `state`: Integer indicating the observations or states of each community.
#'     The sequence of states of a given community forms a trajectory.
#'     + `sp1, ..., sp12`: Vectors containing species abundances for each community
#'     state.
#' * `state_dissim`: Object of class `dist` containing Bray-Curtis dissimilarities
#' between every pair of states in `abundance` (see Details).
#' * `segment_dissim`: Object of class `dist` containing the dissimilarities
#' between every pair of trajectory segments in `abundance` (see Details).
#' * `traj_dissim`: Object of class `dist` containing the dissimilarities
#' between every pair of community trajectories in `abundance` (see Details).
#'
#' @details
#' Artificial data was generated following the procedure explained in Box 1 in
#' Sánchez-Pinillos et al. In particular, the initial state of each community was
#' defined using a hypothetical environmental space with optimal locations for
#' 12 species ([coenocliner::coenocline()]). Community dynamics were simulated
#' using a general Lotka-Volterra model ([primer::lvcompg()]).
#'
#' State dissimilarities were calculated through [vegan::vegdist()]. Segment and
#' trajectory dissimilarities were calculated through [ecotraj::segmentDistances()]
#' and [ecotraj::trajectoryDistances()], respectively
#'
#'
"EDR_data"




