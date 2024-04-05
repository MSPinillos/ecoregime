#' Position of a state with respect to a trajectory
#'
#' @description
#' Define the position of a state with respect to a reference trajectory based on
#' its distance from the trajectory and the length and direction of the trajectory.
#'
#' @param d Either a symmetric matrix or an object of class [`dist`] containing
#' the dissimilarities between each pair of states.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `d` belongs.
#' @param states Vector of integers indicating the order of the states in `d` for
#' each trajectory (assign 1 if the state does not belong to any trajectory).
#' @param target_states Vector of integers indicating the indices in `trajectories`
#' and `states` of the ecological states for which their relative position will
#' be calculated.
#' @param reference Vector of the same class of `trajectories` or object of class
#' `RETRA` indicating the reference trajectory to calculate the relative position
#' of the `target_states`
#' @param method Method to calculate the distance and relative position of the
#' `target_states` and the `reference`. One of `"nearest_state"`, `"projection"`
#' or `"mixed"` (see Details).
#' @param coordStates Matrix containing the coordinates of each state (rows) and
#' axis (columns) of a metric ordination space (see Details)
#'
#' @return
#' The function `state_to_trajectory()` returns a data frame of four columns
#' including the `distance` and `relative_position` between the `target_state` and
#' the `reference`.
#' * Depending on the `method`, `distance` is calculated as the dissimilarity
#' between the `target_states` and their respective nearest state in the `reference`
#' or the dissimilarity to their projections onto the `reference`.
#' * The `relative_position` is a value that ranges between 0 (if the nearest
#' state or projected point coincides with the first `reference` state) and 1
#' (if the nearest state or projected point coincides with the last `reference`
#' state).
#'
#' @details
#' `state_to_trajectory()` can calculate the distance and relative position of
#' one or more `target_states` relative to a `reference` trajectory by three
#' different methods:
#' * `"nearest_state"` returns the dissimilarity of the `target_states` to the
#' nearest state of the `reference` trajectory (`distance`) and calculates the
#' relative position of the nearest state within the `reference`.
#' * `"projection"` returns the dissimilarity of the `target_states` to their
#' projection onto the `reference` trajectory and calculates the relative position
#' of the projected state within the `reference`. This method requires `d` to be
#' metric (i.e. to satisfy the triangle inequality). If `d` is not metric,
#' `state_to_trajectory()` calculates the Euclidean distance within a transformed
#' space generated through multidimensional scaling (Borg and Groenen, 2005). To
#' use the state coordinates in a different metric space, use the `coordStates`
#' argument. When the `target_states` cannot be projected onto any of the segments
#' forming the `reference` trajectory, `state_to_trajectory()` returns `NA` for
#' both `distance` and `relative_position`.
#' * `"mixed"` calculates the dissimilarity between the `target_states` and the
#' `reference` trajectory, as well as their relative position by computing its
#' projection onto any of the segments of the reference (analogous to
#' `method = "projection"`). For the `target_states` that cannot be projected,
#' `state_to_trajectory()` uses the nearest state in the `reference` to compute
#' `distance` and `relative_position` (analogous to `method = "nearest_state"`).
#'
#' @author Martina Sánchez-Pinillos
#'
#' @export
#'
#' @examples
#' # State dissimilarities
#' d <- vegan::vegdist(EDR_data$EDR3$abundance[, paste0("sp", 1:12)], method = "bray")
#' trajectories <- EDR_data$EDR3$abundance$traj
#' states <- EDR_data$EDR3$abundance$state
#'
#' # Calculate the representative trajectories of an EDR to be used as reference
#' RT <- retra_edr(d,
#'                trajectories = trajectories,
#'                states = states,
#'                minSegs = 10)
#'
#' # Define the target states
#' target_states <- as.integer(c(1, 16, 55))
#'
#' # Calculate the position of the target states with respect to  the representative
#' # trajectories of an EDR
#' state_to_trajectory(d, trajectories = trajectories,
#'                     states = states,
#'                     target_states = target_states,
#'                     reference = RT,
#'                     method = "nearest_state")
#'
state_to_trajectory <- function(d, trajectories, states, target_states, reference,
                                method, coordStates = NULL) {

  # due to NSE notes in R CMD check
  d1 = d2 = dseg = P = A = st_to_seg = segment = d1_dseg = d2_dseg = relative_position = NULL

  ## WARNING MESSAGES ----------------------------------------------------------

  method <- match.arg(method, c("nearest_state", "projection", "mixed"))

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !inherits(d, "dist")) |
      nrow(as.matrix(d)) != ncol(as.matrix(d))) {
    stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  }

  if (length(trajectories) != nrow(as.matrix(d))) {
    stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
  }
  if (inherits(reference, "RETRA")) {
    if (any(length(grep("-", trajectories)) > 0,
            length(grep("\\[", trajectories)) > 0,
            length(grep("\\]", trajectories)) > 0)) {
      stop("Avoid using '-', '[', ']' in the values of 'trajectories'.")
    }
  }

  if (length(states) != nrow(as.matrix(d))) {
    stop("The length of 'states' must be equal to both dimensions in 'd'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check the format of target_states
  if (!is.integer(target_states) | any(target_states > length(trajectories))) {
    stop("'target_states' needs to be a vector of class 'integer' with values lower than the length of 'trajectories'.")
  }

  # Check the format of reference
  if (!inherits(reference, "RETRA")) {
    if (!all(reference %in% trajectories)) {
      stop("'reference' must be included in 'trajectories'.")
    }
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

  ## REFERENCE STATES ----------------------------------------------------------

  # reference is an object of class "RETRA"
  if (inherits(reference, "RETRA")) {

    # Representative segments
    RT_segments <- lapply(reference, "[", "Segments")
    RT_states <- lapply(RT_segments, function(segs){
      seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segs[[1]])), "-")
      unlist(lapply(seg_components, function(iseg){
        c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
      }))
    })

    # Check that the states of RT are included in d
    RT_in_d <- sapply(RT_states, function(rt){
      all(rt %in% traj_st)
    })
    if (all(RT_in_d == F)) {
      stop("All states in 'reference' must be included in 'd' and specified in 'trajectories' and 'states'.")
    }

    # Remove duplicated states if necessary
    RT_states <- lapply(RT_states, function(iRT){
      iRT[sapply(1:(length(iRT)-1), function(iseg){
        iRT[iseg] != iRT[(iseg+1)]
      })]
    })

    # Indices of the states forming representative trajectories
    ref_states <- lapply(RT_states, function(iRT){
      match(iRT, traj_st)
    })

  }

  # reference is NOT an object of class "RETRA"
  if (!inherits(reference, "RETRA")) {
    ref_states <- lapply(setNames(reference, reference), function(iRT){
      ref <- which(trajectories %in% iRT)
      ref[order(states[ref])]
    })
  }

  # Dissimilarities between consecutive states in reference
  d_ref <- lapply(ref_states, function(iRT){
    d_ref_st <- d_mat[iRT, iRT]
    c(0, sapply(1:(length(iRT)-1), function(iseg){
      d_ref_st[iseg, iseg+1]
    }))
  })

  ## NEAREST_STATE -------------------------------------------------------------

  if (method == "nearest_state") {
    # Calculate the dissimilarity to the nearest state and the relative position
    # of that state
    st_to_traj.ls <- lapply(target_states, function(itarget){
      lapply(names(ref_states), function(iRT){
        df <- data.frame(target_state = itarget,
                         reference = iRT,
                         distance = min(d_mat[itarget, ref_states[[iRT]]]))
        which_state <- which.min(d_mat[itarget, ref_states[[iRT]]])
        df$relative_position = sum(d_ref[[iRT]][1:which_state]) / sum(d_ref[[iRT]])
        return(df)
      })
    })

    # Convert into a data frame
    st_to_traj <- do.call(rbind, lapply(st_to_traj.ls, function(itarget){
      do.call(rbind, itarget)
    }))
  }

  ## PROJECTION ----------------------------------------------------------------

  if (method %in% c("projection", "mixed")) {

    # TRIANGLE INEQUALITY ------------------------------------------------------

    # Dissimilarities to check triangle inequality
    d_tar_ref <- lapply(setNames(target_states, target_states), function(itarget){
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

    ## PROJECTION BASED ON DISSIMILARITIES -------------------------------------

    if (is_metric == T) {
      # Position of target states with respect to their projection onto the reference
      st_to_traj.ls <- lapply(seq_along(d_tar_ref), function(itarget){
        lapply(setNames(names(d_tar_ref[[itarget]]), names(d_tar_ref[[itarget]])), function(iRT){
          # Dissimilarity to the projected state
          d_tar_ref[[itarget]][[iRT]][, P := rowSums(d_tar_ref[[itarget]][[iRT]][, c('d1', 'd2', 'dseg')])]
          d_tar_ref[[itarget]][[iRT]][, A := sqrt(P/2 * (P/2 - d1) * (P/2 - d2) * (P/2 - dseg))]
          d_tar_ref[[itarget]][[iRT]][, st_to_seg := 2 * A / dseg]
          d_tar_ref[[itarget]][[iRT]][, segment := 1:.N]
          # Components of d1 and d2 onto dseg (use 'abs' to avoid negative values for lack of precision)
          d_tar_ref[[itarget]][[iRT]][, d1_dseg := sqrt(abs(d1^2 - st_to_seg^2))]
          d_tar_ref[[itarget]][[iRT]][, d2_dseg := sqrt(abs(d2^2 - st_to_seg^2))]
          # Select projections within segment limits and the minimum state-to-projection dissimilarity
          d_tar_ref[[itarget]][[iRT]] <- d_tar_ref[[itarget]][[iRT]][dseg - d1_dseg >= 0 & dseg - d2_dseg >= 0][which.min(st_to_seg)]

          # Relative position
          if (nrow(d_tar_ref[[itarget]][[iRT]]) > 0) {
            d_tar_ref[[itarget]][[iRT]][, relative_position := (sum(d_ref[[iRT]][1:segment]) + d1_dseg) / sum(d_ref[[iRT]])]
          } else {
            d_tar_ref[[itarget]][[iRT]] <- NULL
          }
          return(d_tar_ref[[itarget]][[iRT]])
        })
      })

      # Convert into a data frame
      st_to_traj <- do.call(rbind, lapply(seq_along(st_to_traj.ls), function(itarget){
        do.call(rbind, lapply(names(st_to_traj.ls[[itarget]]), function(iRT){
          if (!is.null(st_to_traj.ls[[itarget]][[iRT]])) {
            df <- data.frame(target_state = target_states[itarget],
                             reference = iRT,
                             distance = st_to_traj.ls[[itarget]][[iRT]]$st_to_seg,
                             relative_position = st_to_traj.ls[[itarget]][[iRT]]$relative_position)
          } else {
            if (method == "projection") {
              df <- data.frame(target_state = target_states[itarget],
                               reference = iRT,
                               distance = NA,
                               relative_position = NA)
            }
            if (method == "mixed") {
              df <- data.frame(target_state = target_states[itarget],
                               reference = iRT,
                               distance = min(d_mat[target_states[itarget], ref_states[[iRT]]]),
                               relative_position = sum(d_ref[[iRT]][1:which.min(d_mat[target_states[itarget], ref_states[[iRT]]])]) / sum(d_ref[[iRT]]))
            }
          }
          return(df)
        }))
      }))
    }

    # PROJECTION IN A TRANSFORMED STATE SPACE ----------------------------------

    if (is_metric == F) {

      # Use a transformed metric space
      if (is.null(coordStates)) {
        warning("The dissimilarity metric used in 'd' does not satisfy the triangle inequality. \nThe state space was transformed using metric multidimensional scaling.")

        # Compute MDS and extract the state coordinates
        mds <- smacof::mds(d_mat, ndim = nrow(d_mat) - 1)
        coordStates <- mds$conf
      } else {
        if (length(trajectories) != nrow(coordStates)) {
          stop("The number of rows of 'coordStates' must be equal to the number of states in 'd'")
        }
      }

      # Calculate the Euclidean distance between the target states and their
      # projections onto the reference trajectory
      st_to_traj.ls <- lapply(setNames(target_states, target_states), function(itarget){
        lapply(setNames(names(ref_states), names(ref_states)), function(iRT){
          # Coordinates and euclidean distances of the reference states
          coordRef <- coordStates[ref_states[[iRT]], ]
          eucl_ref <- as.matrix(dist(coordRef))
          eucl_ref <- c(0, sapply(1:(length(ref_states[[iRT]]) - 1), function(iseg){
            eucl_ref[iseg, iseg+1]
          }))

          # Calculate the projection of the states onto the reference segments
          st_to_seg <- do.call(rbind, lapply(1:(length(ref_states[[iRT]]) - 1), function(iseg){

            # Projection
            v_seg <- coordRef[(iseg + 1), ] - coordRef[iseg, ]
            u_st <- coordStates[itarget, ] - coordRef[iseg, ]
            compv_u <- (u_st %*% v_seg) / norm(v_seg, type = "2")
            st_to_seg <- sqrt(norm(u_st, type = "2")^2 - compv_u^2)

            # Relative position
            relative_position <- (sum(eucl_ref[1:iseg]) + compv_u)/sum(eucl_ref)

            data.frame(target_state = itarget,
                       reference = iRT,
                       distance = st_to_seg,
                       compv_u = compv_u,
                       segment = iseg,
                       eucl_seg = eucl_ref[iseg + 1],
                       relative_position = relative_position)
          }))

          # Select projections within segment limits and the minimum state-to-projection dissimilarity
          st_to_seg <- st_to_seg[which(st_to_seg$compv_u >= 0 &
                                         st_to_seg$compv_u <= st_to_seg$eucl_seg),
                                 c("target_state", "reference", "distance", "relative_position")]
          st_to_seg <- st_to_seg[which.min(st_to_seg$distance), ]

          if (nrow(st_to_seg) == 0){
            if (method == "projection") {
              st_to_seg <- data.frame(target_state = itarget,
                                      reference = iRT,
                                      distance = NA,
                                      relative_position = NA)
            }
            if (method == "mixed") {
              d_euc <- as.matrix(dist(coordStates))
              st_to_seg <- data.frame(target_state = itarget,
                                      reference = iRT,
                                      distance = min(d_euc[itarget, ref_states[[iRT]]]),
                                      relative_position = sum(eucl_ref[1:which.min(d_euc[itarget, ref_states[[iRT]]])]) / sum(eucl_ref))
            }

          }

          return(st_to_seg)
        })
      })

      # Convert into a data frame
      st_to_traj <- do.call(rbind, lapply(st_to_traj.ls, function(itarget){
        do.call(rbind, itarget)
      }))
      row.names(st_to_traj) <- 1:nrow(st_to_traj)
    }
  }

  return(st_to_traj)

}
