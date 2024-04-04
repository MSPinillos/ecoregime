#' Metrics of trajectory deviation with respect to a reference trajectory
#'
#' @description
#' Set of metrics to analyze the deviation of disturbed trajectories from an
#' ecological dynamic regime (EDR) considering a representative trajectory as the
#' reference. These metrics include the resistance to the disturbance, amplitude,
#' recovery, and net change.
#'
#' @name deviation_metrics
#' @aliases resistance amplitude recovery net_change
#'
#' @param d Either a symmetric matrix or an object of class [`dist`] containing the
#' dissimilarities between each pair of states.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `d` belongs.
#' @param states Vector of integers indicating the order of the states in `d` for
#' each trajectory.
#' @param disturbed_trajectories Vector of the same class as `trajectories` indicating
#' the identifier of the disturbed trajectories.
#' @param disturbed_states Vector of integers included in `states`indicating the
#' first state after the release of the disturbance for each value in
#' `disturbed_trajectories`.
#' @param predisturbed_states Vector of integers included in `states` indicating
#' the last undisturbed state of each `disturbed_trajectories`. The previous states
#' to `disturbed_states` are considered by default.
#' @param reference Object of class `RETRA` indicating the representative trajectory
#' taken as the reference to compute the amplitude, recovery, and net_change of
#' the disturbed trajectories (see Details).
#' @param index Method to calculate amplitude, recovery, or net change (`"absolute"`,
#' `"relative"`; see Details).
#' @param method Method to calculate the distance between the `disturbed_states`
#' or `predisturbed_states` and the `reference` trajectory. One of `"nearest_state"`,
#' `"projection"` or `"mixed"` (see Details).
#'
#' @return
#' * `resistance()` returns a data frame of two columns indicating the resistance
#' value (`Rt`) for each `disturbed_trajectory`.
#' * `amplitude()` returns a data frame of three columns indicating the amplitude
#' value (`A_abs`; `A_rel`) for each `disturbed_trajectory` and `reference`.
#' If `index = c("absolute", "relative")`, both values are included in a data
#' frame of four columns.
#' * `recovery()` returns a data frame of four columns indicating the recovery
#' value (`Rc_abs`; `Rc_rel`) for each `disturbed_trajectory`, post-disturbance
#' state (`state`) and `reference`. If `index = c("absolute", "relative")`, both
#' values are included in a data frame of five columns.
#' * `net_change` returns a data frame of four columns indicating the net change
#' value (`NC_abs`; `NC_rel`) for each `disturbed_trajectory`, post-disturbance
#' state (`state`), and `reference`. If `index = c("absolute", "relative")`, both
#' values are included in a data frame of five columns.
#'
#' @details
#'
#' \strong{Resistance (`resistance()`)}
#'
#' *Resistance* captures the immediate impact of the disturbance as a function
#' of the changes in the state variables (Sánchez-Pinillos et al., 2019).
#'
#' \eqn{
#' Rt = 1 - d_{pre,dist}
#' }
#'
#' \strong{Amplitude (`amplitude()`)}
#'
#' *Amplitude* indicates the direction in which the system is displaced during the
#' disturbance in relation to the `reference` (Sánchez-Pinillos et al., 2024).
#' Positive values indicate that the disturbance displaces the system towards the
#' boundaries of the dynamic regime. Negative values indicate that the disturbance
#' displaces the system towards the representative trajectory.
#'
#' Two indices can be calculated:
#'
#' If `index = "absolute"`,
#'
#' \eqn{
#' A = d_{dist,RT} - d_{pre,RT}
#' }
#'
#' If `index = "relative"`,
#'
#' \eqn{
#' A = \frac{d_{dist,RT} - d_{pre,RT}}{d_{pre,dist}}
#' }
#'
#' \strong{Recovery (`recovery()`)}
#'
#' *Recovery* quantifies the ability of the system to evolve towards the `reference`
#' following the relief of the disturbance (if positive) or move in the direction
#' of the boundaries of the dynamic regime (if negative) (Sánchez-Pinillos et al.,
#' 2024).
#'
#' Two indices can be calculated:
#'
#' If `index = "absolute"`,
#'
#' \eqn{
#' Rc = d_{dist,RT} - d_{post,RT}
#' }
#'
#' If `index = "relative"`,
#'
#' \eqn{
#' Rc = \frac{d_{dist,RT} - d_{post,RT}}{d_{dist,post}}
#' }
#'
#' \strong{Net change (`net_change()`)}
#'
#' *Net change* quantifies the proximity of the system to the `reference` relative to
#' the pre-disturbed state (Sánchez-Pinillos et al., 2024). Positive values indicate
#' that the system eventually evolves towards the boundaries of the dynamic regime.
#' Negative values indicate that the system eventually evolves towards the
#' `reference`.
#'
#' Two indices can be calculated:
#'
#' If `index = "absolute"`,
#'
#' \eqn{
#' NC = d_{post,RT} - d_{pre,RT}
#' }
#'
#' If `index = "relative"`,
#'
#' \eqn{
#' NC = \frac{d_{post,RT} - d_{pre,RT}}{d_{pre,post}}
#' }
#'
#' In all cases:
#' * \eqn{d_{pre,RT}} is the dissimilarity between the `predisturbed_states` and
#' the `reference`.
#' * \eqn{d_{dist,RT}} is the dissimilarity between the `disturbed_states` and
#' the `reference`.
#' * \eqn{d_{post,RT}} is the dissimilarity between the states after `disturbed_states`
#' and the `reference`.
#' * \eqn{d_{pre,dist}} is the dissimilarity contained in `d` between the
#' `predisturbed_states` and the `disturbed_states`.
#' * \eqn{d_{dist,post}} is the dissimilarity contained in `d` between the
#' `disturbed_states` and the post-disturbed states.
#' * \eqn{d_{pre,post}} is the dissimilarity contained in `d` between the
#' `predisturbed_states` and the post-disturbed states.
#'
#' \eqn{d_{pre,RT}}, \eqn{d_{dist,RT}}, and \eqn{d_{post,RT}} are calculated using
#' the function [`state_to_trajectory()`] by three different methods:
#'
#' * If `method = "nearest_state"`, \eqn{d_{pre,RT}}, \eqn{d_{dist,RT}}, and
#' \eqn{d_{post,RT}} are calculated as the dissimilarity between the pre-disturbance,
#' disturbed, or post-disturbance states and their nearest state in the  `reference`.
#'
#' * If `method = "projection"`, \eqn{d_{pre,RT}}, \eqn{d_{dist,RT}}, and
#' \eqn{d_{post,RT}} are calculated as the dissimilarity between the pre-disturbance,
#' disturbed, or post-disturbance states and their projection onto the `reference`.
#'
#' * If `method = "mixed"`, \eqn{d_{pre,RT}}, \eqn{d_{dist,RT}}, and \eqn{d_{post,RT}}
#' are calculated in the same way than `method = "projection"` whenever the
#' pre-disturbance, disturbed and post-disturbance states can be projected onto
#' any segment of the `reference`. Otherwise, \eqn{d_{pre,RT}}, \eqn{d_{dist,RT}},
#' and \eqn{d_{post,RT}} are calculated using the nearest state of the `reference`.
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' Sánchez-Pinillos, M., Leduc, A., Ameztegui, A., Kneeshaw, D., Lloret, F., & Coll, L.
#' (2019). Resistance, resilience or change: Post-disturbance dynamics of boreal
#' forests after insect outbreaks. *Ecosystems* 22, 1886-1901
#' https://doi.org/10.1007/s10021-019-00378-6
#'
#' Sánchez-Pinillos, M., Dakos, V., & Kéfi, S. (2024). Ecological dynamic regimes:
#' A key concept for assessing ecological resilience. *Biological Conservation*
#' 289, 110409 https://doi.org/10.1016/j.biocon.2023.110409
#'
#' @seealso
#' [`retra_edr()`] to identify representative trajectories in an ecological dynamic
#' regime.
#'
#' [`define_retra()`] to generate an object of class`RETRA`.
#'
#' [`state_to_trajectory()`] to calculate the position of a state with respect to
#' a trajectory.
#'
#' @export
#'
#' @examples
#' # Identify the representative trajectories of the EDR from undisturbed trajectories
#' RT <- retra_edr(d = EDR_data$EDR3$state_dissim,
#'                 trajectories = EDR_data$EDR3$abundance$traj,
#'                 states = as.integer(EDR_data$EDR3$abundance$state),
#'                 minSegs = 5)
#'
#' # Abundance matrix including disturbed and undisturbed trajectories
#' abundance <- rbind(EDR_data$EDR3$abundance,
#'                    EDR_data$EDR3_disturbed$abundance, fill = TRUE)
#'
#' # State dissimilarities (Bray-Curtis) for disturbed and undisturbed trajectories
#' d <- vegan::vegdist(abundance[, paste0("sp", 1:12)], method = "bray")
#'
#' # Resistance
#' Rt <- resistance(d = d, trajectories = abundance$traj, states = abundance$state,
#'                  disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
#'                  disturbed_states = abundance[disturbed_states == 1]$state)
#'
#' # Amplitude
#' A <- amplitude(d = d, trajectories = abundance$traj, states = abundance$state,
#'                disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
#'                disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)
#'
#' # Recovery
#' Rc <- recovery(d = d, trajectories = abundance$traj, states = abundance$state,
#'                disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
#'                disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)
#'
#' # Net change
#' NC <- net_change(d = d, trajectories = abundance$traj, states = abundance$state,
#'                  disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
#'                  disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)
#'

#### RESISTANCE ####

resistance <- function (d, trajectories, states, disturbed_trajectories, disturbed_states,
                        predisturbed_states = disturbed_states - 1) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }
  if (length(states) != nrow(as.matrix(d))) {
    stop("The length of 'states' must be equal to both dimensions in 'd'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check the format of disturbed_trajectories, disturbed_states, and predisturbed_states
  if (length(disturbed_trajectories) != length(unique(disturbed_trajectories))) {
    stop("Each value in 'disturbed_trajectories' must be specified once.")
  }
  if (!all(disturbed_trajectories %in% trajectories)) {
    stop("All 'disturbed_trajectories' must be included in 'trajectories'.")
  }
  if (!all(paste0(disturbed_trajectories, "_", disturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'disturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (!all(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'predisturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (any(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
          paste0(disturbed_trajectories, "_", disturbed_states))) {
    stop("A state cannot be included in both 'predisturbed_states' and 'disturbed_states' if it refers to the same trajectory.")
  }

  ## TRAJECTORY-STATE ----------------------------------------------------------

  # Trajectory-state
  traj_st <- paste0(trajectories, "_", states)

  # Convert d into matrix
  if (!is.matrix(d)) {
    d_mat <- as.matrix(d)
  } else {
    d_mat <- d
  }

  ## DISTURBED TRAJECTORIES ----------------------------------------------------

  # Indices of the pre-disturbance and disturbed states
  idist <- match(paste0(disturbed_trajectories, "_", disturbed_states), traj_st)
  ipre <- match(paste0(disturbed_trajectories, "_", predisturbed_states), traj_st)

  ## RESISTANCE ----------------------------------------------------------------

  # Resistance
  Rt <- data.frame(disturbed_trajectories = trajectories[ipre],
                   Rt = sapply(seq_along(ipre), function(idisttraj){
                     1 - d_mat[ipre[[idisttraj]], idist[[idisttraj]]]
                   }))

  return(Rt)

}

################################################################################

#' @rdname deviation_metrics
#' @export

#### AMPLITUDE ####

amplitude <- function (d, trajectories, states, disturbed_trajectories, disturbed_states,
                       predisturbed_states = disturbed_states - 1, reference,
                       index = c("absolute", "relative"), method = "nearest_state") {

  # due to NSE notes in R CMD check
  d1 = d2 = d_dist_RT = d_pre_RT = distance = dseg = target_state = NULL

  ## WARNING MESSAGES ----------------------------------------------------------

  index <- match.arg(index, c("absolute", "relative"), several.ok = T)
  method <- match.arg(method, c("nearest_state", "projection", "mixed"))

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }
  if (any(length(grep("-", trajectories)) > 0,
          length(grep("\\[", trajectories)) > 0,
          length(grep("\\]", trajectories)) > 0)) {
    stop("Avoid using '-', '[', ']' in the values of 'trajectories'.")
  }
  if (length(states) != nrow(as.matrix(d))) {
    stop("The length of 'states' must be equal to both dimensions in 'd'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check the format of disturbed_trajectories, disturbed_states, and predisturbed_states
  if (length(disturbed_trajectories) != length(unique(disturbed_trajectories))) {
    stop("Each value in 'disturbed_trajectories' must be specified once.")
  }
  if (!all(disturbed_trajectories %in% trajectories)) {
    stop("All 'disturbed_trajectories' must be included in 'trajectories'.")
  }
  if (!all(paste0(disturbed_trajectories, "_", disturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'disturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (!all(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'predisturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (any(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
          paste0(disturbed_trajectories, "_", disturbed_states))) {
    stop("A state cannot be included in both 'predisturbed_states' and 'disturbed_states' if it refers to the same trajectory.")
  }

  # Check the format of reference
  if (!inherits(reference, "RETRA")) {
    stop("'reference' must be an object of class 'RETRA'.")
  }

  ## TRAJECTORY-STATE ----------------------------------------------------------

  # Trajectory-state
  traj_st <- paste0(trajectories, "_", states)

  # Convert d into matrix
  if (!is.matrix(d)) {
    d_mat <- as.matrix(d)
  } else {
    d_mat <- d
  }

  ## REPRESENTATIVE TRAJECTORIES -----------------------------------------------

  # Representative segments
  RT_segments <- lapply(reference, "[", "Segments")
  RT_states <- lapply(RT_segments, function(segs){
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segs[[1]])), "-")
    unlist(lapply(seg_components, function(iseg){
      c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
    }))
  })

  # Check that the states of all RT are included in d
  RT_in_d <- sapply(RT_states, function(rt){
    all(rt %in% traj_st)
  })
  if (any(RT_in_d == F)) {
    stop("All states in 'reference' must be included in 'd' and specified in 'trajectories' and 'states'.")
  }

  # Indices of the states forming representative trajectories
  ref_states <- lapply(RT_states, function(iRT){
    match(iRT, traj_st)
  })

  ## DISTURBED TRAJECTORIES ----------------------------------------------------

  # Indices of the pre-disturbance and disturbed states
  idist <- match(paste0(disturbed_trajectories, "_", disturbed_states), traj_st)
  ipre <- match(paste0(disturbed_trajectories, "_", predisturbed_states), traj_st)

  ## STATE DISSIMILARITIES -----------------------------------------------------

  # If method = "projection" and d does not fit triangle inequality, extract
  # state coordinates in MDS
  if (method %in% c("projection", "mixed")) {

    # Check triangle inequality
    all_target <- unique(unlist(c(ipre, idist)))
    d_ref <- lapply(ref_states, function(iRT){
      d_ref_st <- d_mat[iRT, iRT]
      c(0, sapply(1:(length(iRT)-1), function(iseg){
        d_ref_st[iseg, iseg+1]
      }))
    })
    # Dissimilarities to check triangle inequality and calculate projections
    d_tar_ref <- lapply(setNames(all_target, all_target), function(itarget){
      lapply(setNames(names(ref_states), names(ref_states)), function(iRT){
        # Dissimilarities between the target state and all reference states
        d_target_ref <- d_mat[itarget, ref_states[[iRT]]]
        # Compile all dissimilarities in a data.table
        dt <- data.table::data.table(d1 = d_target_ref[1:(length(d_target_ref)-1)],
                                     d2 = d_target_ref[2:length(d_target_ref)],
                                     dseg = d_ref[[iRT]][-1])
        # Check the triangle inequality
        dt[, is_metric := ifelse(d1 + d2 >= dseg & d1 + dseg >= d2 & d2 + dseg >= d1, T, F)]
      })
    })

    # Check if the triangle inequality is satisfied
    is_metric <- all(data.table::rbindlist(lapply(d_tar_ref, data.table::rbindlist))$is_metric == T)

    if (is_metric == F) {
      warning("The dissimilarity metric used in 'd' does not satisfy the triangle inequality. \nThe state space was transformed using metric multidimensional scaling.")
      # Compute MDS and extract the state coordinates
      mds <- smacof::mds(d_mat, ndim = nrow(d_mat) - 1)
      coordStates <- mds$conf
      d_euc <- as.matrix(dist(coordStates))
    }
  }

  # pre_dist, pre_RT, dist_RT
  if (method %in% c("projection", "mixed") && is_metric == F) {
    #pre_dis
    pre_dist <- data.table::data.table(disturbed_trajectories = trajectories[ipre],
                                       d_pre_dist = sapply(seq_along(ipre), function(idisttraj) {
                                         d_euc[ipre[idisttraj], idist[idisttraj]]
                                       }))
    #pre_RT
    pre_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                         target_states = ipre, reference = reference,
                                                         method = method))
    pre_RT[, disturbed_trajectories := trajectories[target_state]]
    pre_RT[, d_pre_RT := distance]

    # dist_RT
    dist_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                          target_states = idist, reference = reference,
                                                          method = method))
    dist_RT[, disturbed_trajectories := trajectories[target_state]]
    dist_RT[, d_dist_RT := distance]

  } else {
    #pre_dist
    pre_dist <- data.table::data.table(disturbed_trajectories = trajectories[ipre],
                                       d_pre_dist = sapply(seq_along(ipre), function(idisttraj) {
                                         d_mat[ipre[idisttraj], idist[idisttraj]]
                                       }))
    # pre_RT
    pre_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                         target_states = ipre, reference = reference,
                                                         method = method))
    pre_RT[, disturbed_trajectories := trajectories[target_state]]
    pre_RT[, d_pre_RT := distance]

    # dist_RT
    dist_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                          target_states = idist, reference = reference,
                                                          method = method))
    dist_RT[, disturbed_trajectories := trajectories[target_state]]
    dist_RT[, d_dist_RT := distance]

  }

  # data.table to include dissimilarity values
  D_values <- merge(merge(pre_RT[, c("disturbed_trajectories", "reference", "d_pre_RT")],
                          dist_RT[, c("disturbed_trajectories", "reference", "d_dist_RT")]),
                    pre_dist, by = "disturbed_trajectories")

  ## COMPUTE THE INDICES -------------------------------------------------------

  # Amplitude
  A_values <- D_values[, c("disturbed_trajectories", "reference")]
  if ("absolute" %in% index) {
    A_values$A_abs <- D_values$d_dist_RT - D_values$d_pre_RT
  }
  if ("relative" %in% index) {
    A_values$A_rel <- (D_values$d_dist_RT - D_values$d_pre_RT)/D_values$d_pre_dist
  }

  return(data.frame(A_values))

}

################################################################################

#' @rdname deviation_metrics
#' @export

#### RECOVERY ####

recovery <- function (d, trajectories, states, disturbed_trajectories, disturbed_states,
                      reference, index = c("absolute", "relative"), method = "nearest_state") {

  # due to NSE notes in R CMD check
  d1 = d2 = d_dist_RT = d_post_RT = distance = dseg = target_state = NULL

  ## WARNING MESSAGES ----------------------------------------------------------

  index <- match.arg(index, c("absolute", "relative"), several.ok = T)
  method <- match.arg(method, c("nearest_state", "projection", "mixed"))

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }
  if (any(length(grep("-", trajectories)) > 0,
          length(grep("\\[", trajectories)) > 0,
          length(grep("\\]", trajectories)) > 0)) {
    stop("Avoid using '-', '[', ']' in the values of 'trajectories'.")
  }
  if (length(states) != nrow(as.matrix(d))) {
    stop("The length of 'states' must be equal to both dimensions in 'd'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check the format of disturbed_trajectories, disturbed_states, and predisturbed_states
  if (length(disturbed_trajectories) != length(unique(disturbed_trajectories))) {
    stop("Each value in 'disturbed_trajectories' must be specified once.")
  }
  if (!all(disturbed_trajectories %in% trajectories)) {
    stop("All 'disturbed_trajectories' must be included in 'trajectories'.")
  }
  if (!all(paste0(disturbed_trajectories, "_", disturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'disturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }

  # Check the format of reference
  if (!inherits(reference, "RETRA")) {
    stop("'reference' must be an object of class 'RETRA'.")
  }

  ## TRAJECTORY-STATE ----------------------------------------------------------

  # Trajectory-state
  traj_st <- paste0(trajectories, "_", states)

  # Convert d into matrix
  if (!is.matrix(d)) {
    d_mat <- as.matrix(d)
  } else {
    d_mat <- d
  }

  ## REPRESENTATIVE TRAJECTORIES -----------------------------------------------

  # Representative segments
  RT_segments <- lapply(reference, "[", "Segments")
  RT_states <- lapply(RT_segments, function(segs){
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segs[[1]])), "-")
    unlist(lapply(seg_components, function(iseg){
      c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
    }))
  })

  # Check that the states of all RT are included in d
  RT_in_d <- sapply(RT_states, function(rt){
    all(rt %in% traj_st)
  })
  if (any(RT_in_d == F)) {
    stop("All states in 'reference' must be included in 'd' and specified in 'trajectories' and 'states'.")
  }

  # Indices of the states forming representative trajectories
  ref_states <- lapply(RT_states, function(iRT){
    match(iRT, traj_st)
  })

  ## DISTURBED TRAJECTORIES ----------------------------------------------------

  # Index of the states in disturbed trajectory
  ids_dist <- lapply(setNames(disturbed_trajectories, disturbed_trajectories), function(idisttraj){
    which(trajectories == idisttraj)
  })

  # Indices of the disturbed and post-disturbance states
  idist <- lapply(setNames(seq_along(ids_dist), names(ids_dist)), function(idisttraj){
    which(traj_st == paste0(disturbed_trajectories[idisttraj], "_", disturbed_states[idisttraj]))
  })

  ipost <- lapply(setNames(seq_along(ids_dist), names(ids_dist)), function(idisttraj){
    all_states <- states[ids_dist[[idisttraj]]]
    i_post <- which(all_states > disturbed_states[[idisttraj]])
    traj_st_post <- paste0(disturbed_trajectories[idisttraj], "_", states[ids_dist[[idisttraj]]][i_post])
    which(traj_st %in% traj_st_post)
  })

  ## STATE DISSIMILARITIES -----------------------------------------------------

  # If method = "projection" and d does not fit triangle inequality, extract
  # state coordinates in MDS
  if (method %in% c("projection", "mixed")) {
    # Check triangle dissimilarity
    all_target <- unique(unlist(c(idist, ipost)))
    d_ref <- lapply(ref_states, function(iRT){
      d_ref_st <- d_mat[iRT, iRT]
      c(0, sapply(1:(length(iRT)-1), function(iseg){
        d_ref_st[iseg, iseg+1]
      }))
    })

    # Dissimilarities to check triangle inequality and calculate projections
    d_tar_ref <- lapply(setNames(all_target, all_target), function(itarget){
      lapply(setNames(names(ref_states), names(ref_states)), function(iRT){
        # Dissimilarities between the target state and all reference states
        d_target_ref <- d_mat[itarget, ref_states[[iRT]]]
        # Compile all dissimilarities in a data.table
        dt <- data.table::data.table(d1 = d_target_ref[1:(length(d_target_ref)-1)],
                                     d2 = d_target_ref[2:length(d_target_ref)],
                                     dseg = d_ref[[iRT]][-1])
        # Check the triangle inequality
        dt[, is_metric := ifelse(d1 + d2 >= dseg & d1 + dseg >= d2 & d2 + dseg >= d1, T, F)]
      })
    })

    # Check if the triangle inequality is satisfied
    is_metric <- all(data.table::rbindlist(lapply(d_tar_ref, data.table::rbindlist))$is_metric == T)

    if (is_metric == F) {
      warning("The dissimilarity metric used in 'd' does not satisfy the triangle inequality. \nThe state space was transformed using metric multidimensional scaling.")
      # Compute MDS and extract the state coordinates
      mds <- smacof::mds(d_mat, ndim = nrow(d_mat) - 1)
      coordStates <- mds$conf
      d_euc <- as.matrix(dist(coordStates))
    }

  }

  # dist_post, dist_RT, post_RT
  if (method %in% c("projection", "mixed") && is_metric == F) {

    # dist_post
    dist_post <- data.table::rbindlist(lapply(seq_along(idist), function(idisttraj) {
      data.table::data.table(disturbed_trajectories = disturbed_trajectories[idisttraj],
                             ipost = ipost[[idisttraj]],
                             d_dist_post = d_euc[idist[[idisttraj]], ipost[[idisttraj]]])
    }))

    # dist_RT
    dist_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                          target_states = unlist(idist), reference = reference,
                                                          method = method))
    dist_RT[, disturbed_trajectories := trajectories[target_state]]
    dist_RT[, d_dist_RT := distance]

    # post_RT
    post_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                          target_states = unlist(ipost), reference = reference,
                                                          method = method))
    post_RT[, disturbed_trajectories := trajectories[target_state]]
    post_RT[, d_post_RT := distance]
    post_RT[, ipost := target_state]


  } else {

    # dist_post
    dist_post <- data.table::rbindlist(lapply(seq_along(idist), function(idisttraj) {
      data.table::data.table(disturbed_trajectories = disturbed_trajectories[idisttraj],
                             ipost = ipost[[idisttraj]],
                             d_dist_post = d_mat[idist[[idisttraj]], ipost[[idisttraj]]])
    }))

    # dist_RT
    dist_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                          target_states = unlist(idist), reference = reference,
                                                          method = method))
    dist_RT[, disturbed_trajectories := trajectories[target_state]]
    dist_RT[, d_dist_RT := distance]

    # post_RT
    post_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                          target_states = unlist(ipost), reference = reference,
                                                          method = method))
    post_RT[, disturbed_trajectories := trajectories[target_state]]
    post_RT[, d_post_RT := distance]
    post_RT[, ipost := target_state]

  }

  # data.table to include dissimilarity values
  D_values <- merge(merge(dist_RT[, c("disturbed_trajectories", "reference", "d_dist_RT")],
                          post_RT[, c("disturbed_trajectories", "reference", "d_post_RT", "ipost")]),
                    dist_post, by = c("disturbed_trajectories", "ipost"))
  D_values[, states := states[ipost]]

  ## COMPUTE THE INDICES -------------------------------------------------------

  # Recovery
  Rc_values <- D_values[, c("disturbed_trajectories", "states", "reference")]
  if ("absolute" %in% index) {
    Rc_values$Rc_abs <- D_values$d_dist_RT - D_values$d_post_RT
  }
  if ("relative" %in% index) {
    Rc_values$Rc_rel <- (D_values$d_dist_RT - D_values$d_post_RT)/D_values$d_dist_post
  }
  data.table::setorderv(Rc_values, c("disturbed_trajectories", "reference", "states"))

  return(data.frame(Rc_values))

}

################################################################################

#' @rdname deviation_metrics
#' @export

#### NET CHANGE ####

net_change <- function (d, trajectories, states, disturbed_trajectories, disturbed_states,
                        predisturbed_states = disturbed_states - 1, reference,
                        index = c("absolute", "relative"), method = "nearest_state") {

  # due to NSE notes in R CMD check
  d1 = d2 = d_post_RT = d_pre_RT = distance = dseg = target_state = NULL

  ## WARNING MESSAGES ----------------------------------------------------------

  index <- match.arg(index, c("absolute", "relative"), several.ok = T)
  method <- match.arg(method, c("nearest_state", "projection", "mixed"))

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }
  if (any(length(grep("-", trajectories)) > 0,
          length(grep("\\[", trajectories)) > 0,
          length(grep("\\]", trajectories)) > 0)) {
    stop("Avoid using '-', '[', ']' in the values of 'trajectories'.")
  }
  if (length(states) != nrow(as.matrix(d))) {
    stop("The length of 'states' must be equal to both dimensions in 'd'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check the format of disturbed_trajectories, disturbed_states, and predisturbed_states
  if (length(disturbed_trajectories) != length(unique(disturbed_trajectories))) {
    stop("Each value in 'disturbed_trajectories' must be specified once.")
  }
  if (!all(disturbed_trajectories %in% trajectories)) {
    stop("All 'disturbed_trajectories' must be included in 'trajectories'.")
  }
  if (!all(paste0(disturbed_trajectories, "_", disturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'disturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (!all(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
           paste0(trajectories, "_", states))) {
    stop("All 'predisturbed_states' of a given trajectory must be included in 'states' for the same trajectory.")
  }
  if (any(paste0(disturbed_trajectories, "_", predisturbed_states) %in%
          paste0(disturbed_trajectories, "_", disturbed_states))) {
    stop("A state cannot be included in both 'predisturbed_states' and 'disturbed_states' if it refers to the same trajectory.")
  }

  # Check the format of reference
  if (!inherits(reference, "RETRA")) {
    stop("'reference' must be an object of class 'RETRA'.")
  }

  ## TRAJECTORY-STATE ----------------------------------------------------------

  # Trajectory-state
  traj_st <- paste0(trajectories, "_", states)

  # Convert d into matrix
  if (!is.matrix(d)) {
    d_mat <- as.matrix(d)
  } else {
    d_mat <- d
  }

  ## REPRESENTATIVE TRAJECTORIES -----------------------------------------------

  # Representative segments
  RT_segments <- lapply(reference, "[", "Segments")
  RT_states <- lapply(RT_segments, function(segs){
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segs[[1]])), "-")
    unlist(lapply(seg_components, function(iseg){
      c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
    }))
  })

  # Check that the states of all RT are included in d
  RT_in_d <- sapply(RT_states, function(rt){
    all(rt %in% traj_st)
  })
  if (any(RT_in_d == F)) {
    stop("All states in 'reference' must be included in 'd' and specified in 'trajectories' and 'states'.")
  }

  # Indices of the states forming representative trajectories
  ref_states <- lapply(RT_states, function(iRT){
    match(iRT, traj_st)
  })

  ## DISTURBED TRAJECTORIES ----------------------------------------------------

  # Index of the states in disturbed trajectory
  ids_dist <- lapply(setNames(disturbed_trajectories, disturbed_trajectories), function(idisttraj){
    which(trajectories == idisttraj)
  })

  ipost <- lapply(setNames(seq_along(ids_dist), names(ids_dist)), function(idisttraj){
    all_states <- states[ids_dist[[idisttraj]]]
    i_post <- which(all_states > disturbed_states[[idisttraj]])
    traj_st_post <- paste0(disturbed_trajectories[idisttraj], "_", states[ids_dist[[idisttraj]]][i_post])
    which(traj_st %in% traj_st_post)
  })

  ipre <- sapply(seq_along(disturbed_trajectories), function(iRT){
    which(traj_st %in% paste0(disturbed_trajectories, "_", predisturbed_states)[iRT])
  })

  ## STATE DISSIMILARITIES -----------------------------------------------------

  # If method = "projection" and d does not fit triangle inequality, extract
  # state coordinates in MDS
  if (method %in% c("projection", "mixed")) {
    # Check triangle dissimilarity
    all_target <- unique(c(ipre, unlist(ipost)))
    d_ref <- lapply(ref_states, function(iRT){
      d_ref_st <- d_mat[iRT, iRT]
      c(0, sapply(1:(length(iRT)-1), function(iseg){
        d_ref_st[iseg, iseg+1]
      }))
    })

    # Dissimilarities to check triangle inequality and calculate projections
    d_tar_ref <- lapply(setNames(all_target, all_target), function(itarget){
      lapply(setNames(names(ref_states), names(ref_states)), function(iRT){
        # Dissimilarities between the target state and all reference states
        d_target_ref <- d_mat[itarget, ref_states[[iRT]]]
        # Compile all dissimilarities in a data.table
        dt <- data.table::data.table(d1 = d_target_ref[1:(length(d_target_ref)-1)],
                                     d2 = d_target_ref[2:length(d_target_ref)],
                                     dseg = d_ref[[iRT]][-1])
        # Check the triangle inequality
        dt[, is_metric := ifelse(d1 + d2 >= dseg & d1 + dseg >= d2 & d2 + dseg >= d1, T, F)]
      })
    })

    # Check if the triangle inequality is satisfied
    is_metric <- all(data.table::rbindlist(lapply(d_tar_ref, data.table::rbindlist))$is_metric == T)

    if (is_metric == F) {
      warning("The dissimilarity metric used in 'd' does not satisfy the triangle inequality. \nThe state space was transformed using metric multidimensional scaling.")
      # Compute MDS and extract the state coordinates
      mds <- smacof::mds(d_mat, ndim = nrow(d_mat) - 1)
      coordStates <- mds$conf
      d_euc <- as.matrix(dist(coordStates))
    }

  }

  # pre_post, pre_RT, post_RT
  if (method %in% c("projection", "mixed") && is_metric == F) {
    # pre_post, d is not metric
    pre_post <- data.table::rbindlist(lapply(seq_along(ipre), function(idisttraj) {
      data.table::data.table(disturbed_trajectories = disturbed_trajectories[idisttraj],
                             ipost = ipost[[idisttraj]],
                             d_pre_post = d_euc[ipre[[idisttraj]], ipost[[idisttraj]]])
    }))

    # pre_RT, d is not metric
    pre_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                         target_states = unlist(ipre), reference = reference,
                                                         method = method))
    pre_RT[, disturbed_trajectories := trajectories[target_state]]
    pre_RT[, d_pre_RT := distance]

    # post_RT, d is not metric
    post_RT <- data.table::data.table(state_to_trajectory(d = d_euc, trajectories = trajectories, states = states,
                                                          target_states = unlist(ipost), reference = reference,
                                                          method = method))
    post_RT[, disturbed_trajectories := trajectories[target_state]]
    post_RT[, d_post_RT := distance]
    post_RT[, ipost := target_state]

  } else {
    # pre_post, d is metric
    pre_post <- data.table::rbindlist(lapply(seq_along(ipre), function(idisttraj) {
      data.table::data.table(disturbed_trajectories = disturbed_trajectories[idisttraj],
                             ipost = ipost[[idisttraj]],
                             d_pre_post = d_mat[ipre[[idisttraj]], ipost[[idisttraj]]])
    }))

    # pre_RT, d is metric
    pre_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                         target_states = unlist(ipre), reference = reference,
                                                         method = method))
    pre_RT[, disturbed_trajectories := trajectories[target_state]]
    pre_RT[, d_pre_RT := distance]

    # post_RT, d is metric
    post_RT <- data.table::data.table(state_to_trajectory(d = d_mat, trajectories = trajectories, states = states,
                                                          target_states = unlist(ipost), reference = reference,
                                                          method = method))
    post_RT[, disturbed_trajectories := trajectories[target_state]]
    post_RT[, d_post_RT := distance]
    post_RT[, ipost := target_state]
  }

  # data.table to include dissimilarity values
  D_values <- merge(merge(pre_RT[, c("disturbed_trajectories", "reference", "d_pre_RT")],
                          post_RT[, c("disturbed_trajectories", "reference", "d_post_RT", "ipost")]),
                    pre_post, by = c("disturbed_trajectories", "ipost"))
  D_values[, states := states[ipost]]

  ## COMPUTE THE INDICES -------------------------------------------------------

  # Net change
  NC_values <- D_values[, c("disturbed_trajectories", "states", "reference")]
  if ("absolute" %in% index) {
    NC_values$NC_abs <- D_values$d_post_RT - D_values$d_pre_RT
  }
  if ("relative" %in% index) {
    NC_values$NC_rel <- (D_values$d_post_RT - D_values$d_pre_RT)/D_values$d_pre_post
  }
  data.table::setorderv(NC_values, c("disturbed_trajectories", "reference", "states"))

  return(data.frame(NC_values))

}
