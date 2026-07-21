#' Predicted Ecological Trajectories in Ecological Dynamic Regimes (PETRA-EDR)
#'
#' @description
#' Predicts ecological trajectories from a target state or a sequence of states
#' based on their k-nearest states in the trajectories forming a reference ecological
#' dynamic regime (EDR) and the previous and following states in their respective
#' trajectories.
#'
#' @param state_var Object containing the state variables for each trajectory
#' state in both the reference dynamic regime and the `targets`. The class must
#' match the one required in `d_function`. It can be a list if that is required
#' in `d_function`.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `state_var` belongs.
#' @param states Vector of integers indicating the order of the states in `state_var`
#' for each trajectory.
#' @param targets Vector indicating the trajectory or site of the target or
#' targets for which longer trajectories must be predicted.
#' @param d_function Either a function or a non-empty character string naming the
#' function to be called to compute state dissimilarities (see [`do.call`]). The
#' output of `d_function` must be an object of class [`dist`].
#' @param d_args A list of arguments to the `d_function` call (see [`do.call`]).
#' @param d Either a symmetric matrix or an object of class [`dist`] containing the
#' dissimilarities between each pair of states in `state_var`. The elements need
#' to follow the order indicated by `trajectories` and `states`. If `NULL`, `d`
#' is calculated by calling `d_function`.
#' @param k Vector of integers indicating the number of near states to the target
#' in the trajectories composing the reference ecological dynamic regime.
#' @param eps Numeric vector indicating the dissimilarity threshold for  each target
#' beyond which trajectories in the dynamic regime are not considered.
#' @param minPts Vector of integers including the minimum number of states required
#' to predict preceding and following states for each target.
#' @param w_function String indicating the weighting function to estimate the state
#' variables of the predicted states for each target: `"linear"`, `"power"`,
#' `"exponential"`, `"Gaussian"`, `"hyperbolic"`, `"spherical"`.
#' @param alpha Numeric value of the shape parameter used by `w_function` (required
#' for `"power"`, `"exponential"`, and `"hyperbolic"`).
#' @param w Numeric vector of the same length as `trajectories` giving the weights
#' to estimate the state variables in the predicted trajectory. If `targets`
#' is a vector of two or more elements, list of the same length than `targets`
#' containing a numeric vector with trajectory weights for each target.
#' @param method String naming the method used to compute the states forming the
#' predicted trajectory: `"mean"` (default) or `"medoid"`.
#' @param direction String or integer indicating whether the prediction is done
#' backward (`"backward"` or `-1`), forward (`"forward"` or `1`), or in both
#' directions (`"both"` or `2`, default).
#' @param return_args Logic value indicating whether a list containing the arguments
#' provided in the function call must be returned.
#'
#' @details
#' The PETRA-EDR algorithm (`petra_edr()`) is based on the search for the `k`
#' nearest states of the trajectories forming a reference EDR to the target and
#' a moving window towards both directions of their corresponding trajectories.
#'
#' First, PETRA-EDR looks for the target's `k` nearest states in the trajectories
#' forming the reference EDR. Alternatively, PETRA-EDR may consider all states in
#' a pre-defined radius `eps` by assigning `k` the total number of states in the
#' trajectories of the reference EDR.
#'
#' Then, consecutive states forming the predicted trajectory are defined based
#' on the previous and following states of the `k` nearest states in their respective
#' trajectories.
#'
#' Predicted states are defined by averaging the values of the state variables
#' of at least `minPts` states in the trajectories to which the k-nearest states
#' belong (`method = "mean"`) or using the state variables of the medoid of those
#' states (`method = "medoid"`).
#'
#' To estimate a weighted mean of the state variables, it is necessary to indicate
#' the weight (`w`) that must be used in each trajectory of the EDR. Different
#' weights can be used for each target by providing a list of numeric vectors
#' with the weights assigned to each trajectory state. Alternatively, one can
#' specify the weighting function to be applied depending on the dissimilarity
#' between the target and the states in the trajectories forming the EDR:
#'
#'
#' **`"linear"`**
#'
#' \deqn{w(d_i) = 1 - \frac{d_i}{d_{max}}}
#'
#' **`"power"`**
#'
#' \deqn{w(d_i) = 1 - {(\frac{d_i}{d_{max}})^{\alpha}}}
#'
#' **`"exponential"`**
#'
#' \deqn{w(d_i) = e^{\frac{- \alpha d_i}{d_{max}}}}
#'
#' **`"Gaussian"`**
#'
#' \deqn{w(d_i) = e^{-(\frac{d_i}{d_{max}})^2}}
#'
#' **`"hyperbolic"`**
#'
#' \deqn{w(d_i) = \frac{(1 + \frac{d_i}{d_{max}})^{- \alpha} - 2^{- \alpha}}{1 - 2^{- \alpha}}}
#'
#' **`"spherical"`**
#'
#' \deqn{w(d_i) = 1 - 1.5 \frac{d_i}{d_{max}} + 0.5(\frac{d_i}{d_{max}})^3}
#'
#' where \eqn{d_i} is the dissimilarity between the target and the trajectory \eqn{i},
#' \eqn{d_{max}} is the maximum dissimilarity between the target and all trajectories,
#' and \eqn{\alpha} is a shape parameter.
#'
#' @returns
#' The function `petra_edr()` returns an object of class `PETRA`, which is a list
#' with elements containing the following information:
#' \describe{
#' \item{`arguments`}{List of arguments as specified when the function is called
#' (only if `return_args = TRUE`).}
#' \item{`k_dist`}{Data frame of five columns including `target`: identifier of
#' the target provided in the argument `targets`; `target_state`: indices
#' of the states in the extremes of the `targets`; `k_trajectories`: vector
#' indicating the trajectory or site to which each of the k-nearest states belongs;
#' `k_states`: vector indicating the order of the k-nearest states in their trajectory
#' (`k_trajectories`); `k_dist`: dissimilarity between the target and each of
#' the k-nearest states.}
#' \item{`predicted_dist`}{Data frame of six columns including, for each predicted
#' state (`predicted_state`) of the target (`target`) the number of states
#' in the EDR used in the averaging process to calculate the predicted states
#' (`N`; \eqn{N \ge MinPts}) and the mean (`mean_dist`), standard deviation
#' (`sd_dist`), minimum (`min_dist`), and maximum (`max_dist`) dissimilarities to
#' the predicted states.}
#' \item{`state_var`}{Updated `state_var` object containing the state variables
#' of the states forming the predicted trajectory/ies.}
#' \item{`trajectories`}{Vector indicating the target to which each
#' state in the predicted `state_var` object belongs.}
#' \item{`states`}{Vector of integers indicating the order of the states in the
#' predicted `state_var` object for each target. In order to keep the original
#' state value of the targets, there can be values lower or equal to zero.}
#' }
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
#' Forecasting ecological trajectories from ecological dynamic regimes to improve
#' resilience analysis. *Methods in Ecology and Evolution* (in press).
#'
#' @seealso
#' [`MPD()`] for estimating the prediction accuracy of `petra-edr()` outputs.
#'
#' [`plot.PETRA()`] for plotting predicted trajectories in an ordination space
#' representing the state space of the EDR.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   # Example 1 -----------------------------------------------------------------
#'   # Compute the predicted trajectory of a target composed of one state
#'
#'   # State variables for the states in the trajectories forming the EDR
#'   EDR_var <- EDR_data$EDR1$abundance
#'
#'   # State variables of the target
#'   target_var <- data.table::data.table(sp1 = 30, sp2 = 10, sp3 = 13, sp4 = 12,
#'                                        sp5 = 5, sp6 = 2, sp7 = 6, sp8 = 3,
#'                                        sp9 = 7, sp10 = 8, sp11 = 2, sp12 = 3)
#'
#'   # Define state_var including the state variables for the states in the EDR and
#'   # the target
#'   state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
#'
#'   # Define the function used to calculate state dissimilarities and its arguments
#'   # For example, the Canberra dissimilarity.
#'   d_function = "vegan::vegdist"
#'   d_args = list(x = state_var, method = "canberra")
#'
#'   # Compute PETRA-EDR
#'   petra <- petra_edr(state_var = state_var,
#'                      trajectories = c(EDR_var$traj, "target"),
#'                      states = as.integer(c(EDR_var$state, 1)),
#'                      targets = "target",
#'                      k = 5L,
#'                      minPts = 2L,
#'                      d_function = d_function,
#'                      d_args = d_args)
#'
#'   # Example 2 -----------------------------------------------------------------
#'   # Compute the predicted trajectory of two targets using different parameter
#'   # values
#'
#'   # State variables for the states in the trajectories forming the EDR
#'   EDR_var <- EDR_data$EDR1$abundance
#'
#'   # State variables of the target states
#'   target_var <- data.table::data.table(
#'     traj = c("target1", "target1", "target2"), state = c(1, 2, 1),
#'     sp1 = c(30, 31, 3), sp2 = c(10, 9, 70), sp3 = c(13, 8, 3), sp4 = c(12, 9, 4),
#'     sp5 = c(5, 4, 4), sp6 = c(2, 3, 3), sp7 = c(6, 7, 2), sp8 = c(3, 5, 2),
#'     sp9 = c(7, 6, 3), sp10 = c(8, 7, 3), sp11 = c(2, 3, 2), sp12 = c(3, 3, 2))
#'
#'   # Define state_var including the state variables for the states in the EDR and
#'   # the targets
#'   state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
#'                                 target_var[, -c("traj", "state")]))
#'
#'   # Define trajectories and states in state_var and the ID of the targets
#'   trajectories <- c(EDR_var$traj, target_var$traj)
#'   states <- as.integer(c(EDR_var$state, target_var$state))
#'   targets <- c("target1", "target2")
#'
#'   # Define the function used to calculate state dissimilarities and its arguments
#'   d_function = "vegan::vegdist"
#'   d_args <- list(x = state_var, method = "bray")
#'
#'   # Compute PETRA-EDR
#'   petra <- petra_edr(state_var = state_var,
#'                      trajectories = trajectories,
#'                      states = states,
#'                      targets = targets,
#'                      k = c(5L, 3L),
#'                      minPts = c(2L, 3L),
#'                      eps = c(NA, 0.6),
#'                      d_function = d_function,
#'                      d_args = d_args,
#'                      method = "mean",
#'                      w_function = c("exponential", NA),
#'                      alpha = c(3, NA))
#' }
#'
#'
petra_edr <- function(state_var, trajectories, states, targets,
                      d_function, d_args, d = NULL,
                      k, eps = NULL, minPts = k,
                      w_function = NULL, alpha = 1, w = NULL,
                      method = "mean", direction = 2, return_args = FALSE) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for state_var, trajectories, states
  if (!any(inherits(state_var, "matrix"),
           inherits(state_var, "data.frame"),
           inherits(state_var, "list"))) {
    stop("'state_var' needs to be of any of these classes: 'matrix', 'data.frame', 'list'")
  }
  if (length(trajectories) != length(states)) {
    stop("The length of 'trajectories' must be equal to the length of 'states'.")
  }
  if (inherits(state_var, "list")) {
    if (length(state_var) != length(trajectories)) {
      stop("The length of 'state_var' must be equal to the length of 'trajectories' and 'states'.")
    }
  } else {
    if (nrow(state_var) != length(trajectories)) {
      stop("The number of elements in 'state_var' must be equal to the length of 'trajectories' and 'states'.")
    }
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check that all targets correspond to an index in trajectories and states
  if (!all(targets %in% trajectories)) {
    stop("'target' needs to be included in 'trajectories.")
  }

  # Length of k, eps, and minPts
  if (length(targets) > 1) {
    if (length(k) == 1) {
      k <- rep(k, length(targets))
    }
    if (length(minPts) == 1) {
      minPts <- rep(minPts, length(targets))
    }
    if (!is.null(eps) & length(eps) == 1) {
      eps <- rep(eps, length(targets))
    }

    if (length(k) < length(targets)) {
      stop(cat("'k' must be a numeric vector of length 1 or", length(targets)))
    }
    if (length(minPts) < length(targets)) {
      stop(cat("'minPts' must be a numeric vector of length 1 or", length(targets)))
    }
    if (!is.null(eps) & length(eps) < length(targets)) {
      stop(cat("'eps' must be a numeric vector of length 1 or", length(targets)))
    }
  }

  # Format of k, eps, and minPts
  if (any(!is.integer(k), !is.integer(minPts))) {
    stop("'k' and 'minPts' need to be integer")
  }
  if (any(minPts > k)) {
    stop("'minPts' cannot be greater than 'k'.")
  }
  if (!is.null(eps) & !is.numeric(eps)) {
    stop("'eps' needs to be numeric.")
  }

  # Check d
  if (!is.null(d)) {
    if (all(!is.matrix(d), !inherits(d, "dist")) |
        nrow(as.matrix(d)) != ncol(as.matrix(d))) {
      stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
    }
  }

  # Check d_function
  if (!all(names(d_args) %in% names(formals(eval(parse(text = d_function)))))) {
    stop(cat("Invalid arguments for 'd_function': ", paste(setdiff(names(d_args), names(formals(d_function))), collapse = ", ")))
  }

  # Check method
  method <- match.arg(method, c("mean", "medoid"))

  # Check weights
  if (method == "mean" & !is.null(w)) {
    if (length(targets) == 1 & is.numeric(w)) {
      if (length(w) != length(trajectories)) {
        stop(cat("'w' must be a numeric vector of length", length(trajectories)))
      }
      if (all(w == 0)) {
        w <- rep(1, length(trajectories))
      }
    } else {
      if (!(is.list(w) & length(w) == length(targets)) & !(is.numeric(w) & length(w) == length(trajectories))) {
        stop(paste("'w' must be a list of length", length(targets), "or a numeric vector of length", length(trajectories)))
      }

      if (is.list(w)) {
        if (length(w) != length(targets)) {
          stop(cat("'w' must be a list of length", length(targets)))
        }
        for (i in seq_along(targets)) {
          if (length(w[[i]]) != length(trajectories) || !is.numeric(w[[i]])) {
            stop(paste("The elements in 'w' must be numeric vectors of lengths:", paste(length(trajectories), collapse = ", ")))
          }
          if (all(w[[i]] == 0)) {
            w[[i]] <- rep(1, length(trajectories))
          }
        }
      } else {
        w <- rep(list(w), length(targets))
      }
    }
  }

  # Check w_function
  if (!is.null(w_function)) {
    if (length(targets) > 1) {
      if (length(w_function) == 1) {
        w_function <- rep(w_function, length(targets))
      }
      if (length(w_function) < length(targets)) {
        stop(cat("'w_function' must be a vector of length 1 or", length(targets)))
      }
    }
    w_function <- sapply(w_function, function(ifun){
      if (is.na(ifun)) {
        NA
      } else {
        match.arg(ifun, c("linear", "power", "exponential", "Gaussian", "hyperbolic", "spherical"))
      }
    })
    if (any(w_function %in% c("power", "exponential", "hyperbolic"))) {
      if (length(alpha) == 1) {
        alpha <- rep(alpha, length(targets))
      }
      if (length(alpha) < length(targets)) {
        stop(cat("'alpha' must be a vector of length 1 or", length(targets)))
      }
    }
  }

  # Check direction
  if (is.numeric(direction) & !direction %in% c(-1, 1, 2)) {
    stop("Wrong input for 'direction'")
  } else if (is.character(direction)){
    direction <- match.arg(direction, c("backward", "forward", "both", -1, 1, 2))
  }

  ## INITIAL SETTINGS ----------------------------------------------------------

  # Give an identifier for each state in the EDR
  ID_states <- paste0(trajectories, "_", states)

  # Identify the argument corresponding to state_var in d_args
  istate_var <- which(vapply(d_args, identical, state_var, FUN.VALUE = logical(1)))

  # Calculate d considering the states in the EDR and the targets
  if (is.null(d)) {
    d <- do.call(eval(parse(text = d_function)), args = d_args)
    if (!inherits(d, "dist")) {
      stop("'d_function' must return an object of class 'dist'.")
    }
  }

  # Dissimilarity matrix
  d_mat <- as.matrix(d)

  # Check that there is the same number of elements in state_var, trajectories,
  # and states
  if (length(trajectories) != nrow(d_mat)) {
    stop("The length of 'trajectories' and 'states' needs to be equal to the number
         of elements in 'state_var' and 'd' (if provided).")
  }

  # Calculate weights
  if (method == "mean" & is.null(w) && !is.null(w_function)) {
    w <- calculate_w(d = d_mat, trajectories = trajectories, states = states, targets = targets,
                     k = k, eps = eps, w_function = w_function, alpha = alpha,
                     w_per_state = T)

    if (is.numeric(w) && all(w == 0)) {
      w <- rep(1, length(w))
    } else if (inherits(w, "list")) {
      for (i in seq_along(w)) {
        if (all(w[[i]] == 0)) {
          w[[i]] <- rep(1, length(w[[i]]))
        }
      }
    }

  }

  # Save the arguments
  arguments = list(trajectories = trajectories, states = as.integer(states), targets = targets,
                   k = k, minPts = minPts, eps = eps,
                   d_function = d_function, d_args = d_args, args_state_var = istate_var,
                   d = as.dist(d), method = method, w = w)

  # Identify different targets
  ID_targets <- which(trajectories %in% targets)

  petra_target <- lapply(setNames(seq_along(unique(targets)), unique(targets)), function(itarget){

    # TARGET SETTINGS ----------------------------------------------------------

    target <- which(trajectories == unique(targets)[itarget])

    # Target data
    if (inherits(state_var, "list")) {
      state_var_target <- state_var[target]
    } else {
      state_var_target <- state_var[target, ]
    }

    trajectories_target <- trajectories[target]
    states_target <- states[target]

    # If there are two or more states in the target, identify the first and last
    if (length(target) > 1) {
      target1 <- target[which.min(states_target)]
      target2 <- target[which.max(states_target)]
    } else {
      target1 <- target
      target2 <- target
    }

    # Define weights if necessary
    if (method == "mean" & !is.null(w)) {
      if (is.numeric(w)) {
        weight <- w
      }
      if (is.list(w)) {
        weight <- w[[itarget]]
      }
    }

    ## BUILD PHASE ---------------------------------------------------------------

    # Find the k nearest states backward and forward
    kns_back <- kns_forw <- NULL
    if (direction %in% c("backward", "both", -1, 2)) {
      kns_back <- na.omit(c(setdiff(order(d_mat[, target1]), ID_targets)[1:k[itarget]], target1))
    }
    if (direction %in% c("forward", "both", 1, 2)) {
      kns_forw <- na.omit(c(setdiff(order(d_mat[, target2]), ID_targets)[1:k[itarget]], target2))
    }

    # k-nearest neighbors in eps
    if (!is.null(eps) && !is.na(eps[itarget])) {
      if (direction %in% c("backward", "both", -1, 2)) {
        kns_back <- kns_back[d_mat[kns_back, target1] < eps[itarget]]
      }
      if (direction %in% c("forward", "both", 1, 2)) {
        kns_forw <- kns_forw[d_mat[kns_forw, target2] < eps[itarget]]
      }

      if (length(kns_back) <= 1 & length(kns_forw) <= 1) {
        warning(cat(paste0("target = ", trajectories[target1], ":"), "There are not states at eps < ", eps[itarget], '\n'))
      }
    }

    # k-nearest neighbors with w > 0
    if (!is.null(w)) {
      if (direction %in% c("backward", "both", -1, 2)) {
        kns_back <- kns_back[weight[kns_back[kns_back != target1]] > 0]
      }
      if (direction %in% c("forward", "both", 1, 2)) {
        kns_forw <- kns_forw[weight[kns_forw[kns_forw != target2]] > 0]
      }

      if (length(kns_back) <= 1 & length(kns_forw) <= 1) {
        warning(cat(paste0("target = ", trajectories[target1], ":"), "There are not states at eps < ", eps, '\n'))
      }
    }


    # Table summarizing dissimilarities with kns
    if (target1 == target2) {
      k_dist <- data.table::data.table(target = trajectories[target1], target_state = states[target1],
                                       k_trajectories = trajectories[unique(c(kns_back, kns_forw))],
                                       k_states = states[unique(c(kns_back, kns_forw))],
                                       k_dist = d_mat[unique(c(kns_back, kns_forw)), target1])[target != k_trajectories, ]

      # The target has no neighbors at d < eps
      if (nrow(k_dist) == 0) {
        k_dist <- data.table::data.table(target = trajectories[target1], target_state = states[target1],
                                         k_trajectories = NA,
                                         k_states = NA,
                                         k_dist = NA)
      }

    } else {
      k_dist <- rbind(data.table::data.table(target = trajectories[target1], target_state = states[target1],
                                             k_trajectories = trajectories[kns_back],
                                             k_states = states[kns_back],
                                             k_dist = d_mat[kns_back, target1]),

                      data.table::data.table(target = trajectories[target2], target_state = states[target2],
                                             k_trajectories = trajectories[kns_forw],
                                             k_states = states[kns_forw],
                                             k_dist = d_mat[kns_forw, target2])
      )[paste0(target, "_", target_state) != paste0(k_trajectories, "_", k_states)]

      # The target has no neighbors at d < eps
      if (nrow(k_dist) == 0) {
        k_dist <- rbind(data.table::data.table(target = trajectories[target1], target_state = states[target1],
                                               k_trajectories = NA,
                                               k_states = NA,
                                               k_dist = NA),

                        data.table::data.table(target = trajectories[target2], target_state = states[target2],
                                               k_trajectories = NA,
                                               k_states = NA,
                                               k_dist = NA))
      }
    }

    # Table identifying the states order
    traj_st <- data.table::data.table(traj = trajectories, st = states)
    data.table::setorder(traj_st, traj, st)
    traj_st[, order := 1:(.N), by = traj]
    traj_st[, ID_states := paste0(traj, "_", st)]
    traj_st[, ID_order := paste0(traj, "_", order)]

    # Identify backward states
    i = 1
    traj_st_back <- traj_st[ID_states %in% paste0(trajectories[kns_back], "_", states[kns_back])]
    ID_back <- traj_st[ID_order %in% paste0(traj_st_back$traj, "_", traj_st_back$order - i)]$ID_states
    iback <- which(ID_states %in% ID_back)

    # Identify forward states
    j = 1
    traj_st_forw <- traj_st[ID_states %in% paste0(trajectories[kns_forw], "_", states[kns_forw])]
    ID_forw <- traj_st[ID_order %in% paste0(traj_st_forw$traj, "_", traj_st_forw$order + j)]$ID_states
    iforw <- which(ID_states %in% ID_forw)

    # List to include the dissimilarities between the predicted states and the trajectory states
    predicted_dist <- lapply(states_target, function(i){
      NA
    })


    ## SWEEP BACKWARD ----------------------------------------------------------

    # Sweep backward
    while (length(iback) >= minPts[itarget]) {

      # Define the next state
      if (method == "mean") {
        if (is.null(w)) {
          if (inherits(state_var, "list")) {
            state_var_new <- Reduce("+", state_var[iback]) / length(state_var[iback])
          } else {
            state_var_new <- state_var[target1, ]
            for (icol in 1:ncol(state_var_new)) {
              state_var_new[, icol] <- mean(state_var[iback, icol])
            }
          }
        } else if (!is.null(w)) {
          if (inherits(state_var, "list")) {
            state_var_new <- lapply(iback, function(ib){
              state_var[[ib]] * weight[ib]
            })
            state_var_new <- Reduce("+", state_var_new) / sum(weight[iback])
          } else {
            state_var_new <- state_var[target1, ]
            for(icol in 1:ncol(state_var_new)) {
              state_var_new[, icol] <- weighted.mean(state_var[iback, icol], weight[iback])
            }
          }
        }
      } else if (method == "medoid") {
        if (length(iback) > 1) {
          imedoid <- which.min(colMeans(d_mat[iback, iback]))
        } else {
          imedoid <- 1
        }
        if (inherits(state_var, "list")) {
          state_var_new <- state_var[[iback[imedoid]]]
        } else {
          state_var_new <- state_var[iback[imedoid], ]
        }
      }

      # Calculate the average distance of each predicted state to the neighboring trajectories
      iback_dist <- numeric()
      if (inherits(state_var, "list")) {
        for (iiback in iback) {
          d_args[[istate_var]] <- list(state_var[[iiback]], state_var_new)
          class(d_args[[istate_var]]) <- class(state_var)
          iiback_dist <- do.call(eval(parse(text = d_function)), args = d_args)
          if (!inherits(iiback_dist, "dist")) {
            stop("'d_function' must return an object of class 'dist'.")
          }
          iback_dist <- c(iback_dist, iiback_dist)
        }

      } else {
        for (iiback in iback) {
          d_args[[istate_var]] <- rbind(state_var[iiback, ], state_var_new)
          iiback_dist <- do.call(eval(parse(text = d_function)), args = d_args)
          if (!inherits(iiback_dist, "dist")) {
            stop("'d_function' must return an object of class 'dist'.")
          }
          iback_dist <- c(iback_dist, iiback_dist)
        }
      }
      predicted_dist <- c(predicted_dist, list(iback_dist))

      # Add the new predicted state to state_var
      if (inherits(state_var, "list")) {
        state_var_target <- c(state_var_target, list(state_var_new))
      } else {
        state_var_target <- rbind(state_var_target, state_var_new)
      }

      trajectories_target <- c(trajectories_target, trajectories[target1])
      states_target <- c(states_target, states[target1] - i)

      i = i + 1

      # Identify backward states
      traj_st_back <- traj_st[ID_states %in% paste0(trajectories[kns_back], "_", states[kns_back])]
      ID_back <- traj_st[ID_order %in% paste0(traj_st_back$traj, "_", traj_st_back$order - i)]$ID_states
      iback <- which(ID_states %in% ID_back)

    }

    ## SWEEP FORWARD -----------------------------------------------------------

    # Sweep forward
    while (length(iforw) >= minPts[itarget]) {

      # Define the next state
      if (method == "mean") {
        if (is.null(w)) {
          if (inherits(state_var, "list")) {
            state_var_new <- Reduce("+", state_var[iforw]) / length(state_var[iforw])
          } else {
            state_var_new <- state_var[target2, ]
            for (icol in 1:ncol(state_var_new)) {
              state_var_new[, icol] <- mean(state_var[iforw, icol])
            }
          }
        } else if (!is.null(w)) {
          if (inherits(state_var, "list")) {
            state_var_new <- lapply(iforw, function(ifw){
              state_var[[ifw]] * weight[ifw]
            })
            state_var_new <- Reduce("+", state_var_new) / sum(weight[iforw])
          } else {
            state_var_new <- state_var[target2, ]
            for (icol in 1:ncol(state_var_new)) {
              state_var_new[, icol] <- weighted.mean(state_var[iforw, icol], weight[iforw])
            }
          }
        }
      } else if (method == "medoid") {
        if (length(iforw) > 1) {
          imedoid <- which.min(colMeans(d_mat[iforw, iforw]))
        } else {
          imedoid <- 1
        }
        if (inherits(state_var, "list")) {
          state_var_new <- state_var[[iforw[imedoid]]]
        } else {
          state_var_new <- state_var[iforw[imedoid], ]
        }
      }

      # Calculate the average distance of each predicted state to the neighbouring trajectories
      iforw_dist <- numeric()
      if (inherits(state_var, "list")) {
        for (iiforw in iforw) {
          d_args[[istate_var]] <- list(state_var[[iiforw]], state_var_new)
          class(d_args[[istate_var]]) <- class(state_var)
          iiforw_dist <- do.call(eval(parse(text = d_function)), args = d_args)
          if (!inherits(iiforw_dist, "dist")) {
            stop("'d_function' must return an object of class 'dist'.")
          }
          iforw_dist <- c(iforw_dist, iiforw_dist)
        }

      } else {
        for (iiforw in iforw) {
          d_args[[istate_var]] <- rbind(state_var[iiforw, ], state_var_new)
          iiforw_dist <- do.call(eval(parse(text = d_function)), args = d_args)
          if (!inherits(iiforw_dist, "dist")) {
            stop("'d_function' must return an object of class 'dist'.")
          }
          iforw_dist <- c(iforw_dist, iiforw_dist)
        }
      }
      predicted_dist <- c(predicted_dist, list(iforw_dist))

      # Add the new predicted state to state_var
      if (inherits(state_var, "list")) {
        state_var_target <- c(state_var_target, list(state_var_new))
      } else {
        state_var_target <- rbind(state_var_target, state_var_new)
      }

      trajectories_target <- c(trajectories_target, trajectories[target2])
      states_target <- c(states_target, states[target2] + j)

      j = j + 1

      # Identify forward states
      traj_st_forw <- traj_st[ID_states %in% paste0(trajectories[kns_forw], "_", states[kns_forw])]
      ID_forw <- traj_st[ID_order %in% paste0(traj_st_forw$traj, "_", traj_st_forw$order + j)]$ID_states
      iforw <- which(ID_states %in% ID_forw)

    }

    ## COMPILE -----------------------------------------------------------------

    # Summarize predicted distances
    predicted_dist <- data.table::rbindlist(lapply(seq_along(predicted_dist), function(ipred){
      data.frame(target = trajectories_target[ipred],
                 predicted_state = states_target[ipred],
                 N = ifelse(is.na(predicted_dist[ipred]), NA, length(predicted_dist[[ipred]])),
                 mean_dist = mean(predicted_dist[[ipred]]),
                 sd_dist = sd(predicted_dist[[ipred]]),
                 min_dist = min(predicted_dist[[ipred]]),
                 max_dist = max(predicted_dist[[ipred]]))
    }))

    return(list(
      state_var_target = state_var_target,
      trajectories_target = trajectories_target,
      states_target = states_target,
      k_dist = k_dist,
      predicted_dist = predicted_dist
    ))

  })

  ## OUTPUT --------------------------------------------------------------------

  # Extract and concatenate state_var_target from petra_target
  if (inherits(state_var, "list")) {
    state_var_target <- do.call("c", lapply(petra_target, "[[", "state_var_target"))
  } else {
    state_var_target <- do.call("rbind", lapply(petra_target, "[[", "state_var_target"))
  }

  trajectories_target <- do.call("c", lapply(petra_target, "[[", "trajectories_target"))
  states_target <- do.call("c", lapply(petra_target, "[[", "states_target"))
  k_dist <- do.call("rbind", lapply(petra_target, "[[", "k_dist"))
  predicted_dist <- do.call("rbind", lapply(petra_target, "[[", "predicted_dist"))

  # Set data in order
  set_order <- order(trajectories_target, states_target)
  if (inherits(state_var, "list")) {
    state_var_target <- state_var_target[set_order]
  } else {
    state_var_target <- state_var_target[set_order, ]
  }
  predicted_dist <- predicted_dist[set_order, ]
  states_target <- states_target[set_order]
  trajectories_target <- trajectories_target[set_order]

  # Keep the original class of state_var
  class(state_var_target) <- class(state_var)

  # Define names based on trajectories and states
  if (inherits(state_var_target, "list")) {
    names(state_var_target) <- paste0(trajectories_target, "_", states_target)
  } else if (inherits(state_var_target, "data.frame")) {
    row.names(state_var_target) <- NULL
  }
  names(trajectories_target) <- NULL
  names(states_target) <- NULL

  # Compile everything in a list
  if (return_args == T) {
    petra <- list(
      arguments = arguments,
      k_dist = k_dist,
      predicted_dist = na.omit(predicted_dist),
      state_var = state_var_target,
      trajectories = trajectories_target,
      states = states_target
    )
  } else {
    petra <- list(
      k_dist = k_dist,
      predicted_dist = na.omit(predicted_dist),
      state_var = state_var_target,
      trajectories = trajectories_target,
      states = states_target
    )
  }

  class(petra) <- "PETRA"

  return(petra)

}

################################################################################

# Declare global variables to ignore them during checks.
utils::globalVariables(c(
  "ID_order", "distance", "k_states", "k_trajectories", "sd", "st", "target", "target_state", "traj"
))

#### CALCULATE WEIGHTS ####

# Auxiliary function to calculate weights as the inverse of the distance
calculate_w <- function(d, trajectories, states, targets, k, eps = NULL,
                        w_function, alpha = NULL, w_per_state = TRUE) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check d
  if (!is.null(d)) {
    if (all(!is.matrix(d), !inherits(d, "dist")) |
        nrow(as.matrix(d)) != ncol(as.matrix(d))) {
      stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
    }
  }

  # Dissimilarity matrix
  d_mat <- as.matrix(d)

  if (length(trajectories) != nrow(d_mat)) {
  stop("The length of 'trajectories' and 'states' needs to be equal to the number
       of elements considered in 'd'.")
  }

  # Check the format for state_var, trajectories, states
  if (length(trajectories) != length(states)) {
    stop("The length of 'trajectories' must be equal to the length of 'states'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check that all targets correspond to an index in trajectories and states
  if (!all(targets %in% trajectories)) {
    stop("'target' needs to be included in 'trajectories.")
  }

  # Format of k, eps, and minPts
  if (!is.integer(k)) {
    stop("'k' need to be integer")
  }
  if (!is.null(eps) & !is.numeric(eps)) {
    stop("'eps' needs to be numeric.")
  }

  # Length of k, eps, and minPts
  if (length(targets) > 1) {
    if (length(k) == 1) {
      k <- rep(k, length(targets))
    }
    if (!is.null(eps) & length(eps) == 1) {
      eps <- rep(eps, length(targets))
    }

    if (length(k) < length(targets)) {
      stop(cat("'k' must be a numeric vector of length 1 or", length(targets)))
    }
    if (!is.null(eps) & length(eps) < length(targets)) {
      stop(cat("'eps' must be a numeric vector of length 1 or", length(targets)))
    }
  }

  # Check w_function
  if (!is.null(w_function)) {
    if (length(targets) > 1) {
      if (length(w_function) == 1) {
        w_function <- rep(w_function, length(targets))
      }
      if (length(w_function) < length(targets)) {
        stop(cat("'w_function' must be a vector of length 1 or ", length(targets)))
      }
    }
    w_function <- sapply(w_function, function(ifun){
      if (is.na(ifun)) {
        NA
      } else {
        match.arg(ifun, c("linear", "power", "exponential", "Gaussian", "hyperbolic", "spherical"))
      }
    })
    if (any(w_function %in% c("power", "exponential", "hyperbolic"))) {
      if (length(alpha) == 1) {
        alpha <- rep(alpha, length(targets))
      }
      if (length(alpha) < length(targets)) {
        stop(cat("'alpha' must be a vector of length 1 or", length(targets)))
      }
    }


  }

  # Identify the targets for which w much be calculated
  id_w <- which(!is.na(w_function))
  targets_w <- targets[id_w]
  k_w <- k[id_w]
  eps_w <- eps[id_w]
  w_function_w <- w_function[id_w]
  alpha_w <- alpha[id_w]

  #-----------------------------------------------------------------------------

  # Identify targets with one or multiple states
  target_length <- vapply(targets_w, function(itarget){
    length(trajectories[trajectories == itarget])
  }, numeric(1))
  target_1 <- targets_w[target_length == 1]
  target_2 <- targets_w[target_length > 1]

  # Identify the states within eps
  if (!is.null(eps)) {
    eps_neighborhood <- lapply(seq_along(targets_w), function(itarget){
      if (targets_w[itarget] %in% target_1) {
        if (!is.na(eps_w[itarget])) {
          unique(trajectories[which(d_mat[, which(trajectories %in% targets_w[itarget])] <= eps_w[itarget])])
        } else {
          unique(trajectories)
        }

      } else if (targets_w[itarget] %in% target_2) {
        if (!is.na(eps_w[itarget])) {
          unique(rep(trajectories,
                     ncol(d_mat[, which(trajectories %in% targets_w[itarget])]))[which(d_mat[, which(trajectories %in% targets_w[itarget])] <= eps_w[itarget])])
        } else {
          unique(trajectories)
        }
      }
    })
  }

  # Identify the k-nearest states
  k_neighbors <- lapply(seq_along(targets_w), function(itarget) {
    if (targets_w[itarget] %in% target_1) {
      # Index of the target
      itarget_1 <- which(trajectories %in% targets_w[itarget])
      neighbors <- unique(trajectories[-itarget_1][order(d_mat[-itarget_1, itarget_1])[1:k_w[itarget]]])
    }

    if (targets_w[itarget] %in% target_2) {
      # Index of the target
      itarget_2 <- which(trajectories %in% targets_w[itarget])
      neighbors <- unique(trajectories[-itarget_2][c(order(d_mat[-itarget_2, itarget_2[1]])[1:k_w[itarget]],
                                                     order(d_mat[-itarget_2, itarget_2[length(itarget_2)]])[1:k_w[itarget]])])
    }

    return(neighbors)
  })

  # Targets with one state -----------------------------------------------------

  if (length(target_1) > 0) {

    # Index of the targets
    itarget_1 <- which(trajectories %in% target_1)

    # Define all trajectories as the reference
    reference <- unique(trajectories[!trajectories %in% target_1])

    # Distance from target_1 to references
    stTraj <- ecoregime::state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                             target_states = itarget_1, reference = reference,
                                             method = "nearest")

    # Target
    stTraj$target <- trajectories[stTraj$target_state]

  }

  # Targets with multiple states -----------------------------------------------

  if (length(target_2) > 0) {

    # Identify trajectories (targets) with one state
    one_st <- names(table(trajectories))[table(trajectories) == 1]

    # Data frame to store dissimilarities
    trTraj <- expand.grid(target = target_2, reference = unique(trajectories[!trajectories %in% one_st]),
                          stringsAsFactors = F)
    trTraj <- trTraj[which(trTraj$target != trTraj$reference), ]

    # Data for the target and each reference trajectory
    sub_trajectories <- lapply(1:nrow(trTraj), function(irow){
      trajectories[trajectories %in% trTraj[irow, ]]
    })
    sub_states <- lapply(1:nrow(trTraj), function(irow){
      states[trajectories %in% trTraj[irow, ]]
    })
    sub_d <- lapply(1:nrow(trTraj), function(irow){
      d_mat[trajectories %in% trTraj[irow, ],
            trajectories %in% trTraj[irow, ]]
    })

    # Calculate the dissimilarity between the target and each trajectory
    trTraj$distance <- vapply(1:nrow(trTraj), function(irow){
      dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(
        d = sub_d[[irow]],
        sites = sub_trajectories[[irow]],
        surveys = sub_states[[irow]]
      ),
      symmetrization = NULL)
      dTraj[row.names(dTraj) %in% trTraj$target[irow], !colnames(dTraj) %in% trTraj$target[irow]]
    }, numeric(1))

  }

  # All targets ----------------------------------------------------------------

  # Compile all targets if necessary
  if (all(length(target_1) > 0, length(target_2) > 0)) {
    w <- rbind(stTraj[, c("target", "reference", "distance")],
               trTraj[, c("target", "reference", "distance")])
  } else if (length(target_1) == 0) {
    w <- trTraj[, c("target", "reference", "distance")]
  } else if (length(target_2) == 0) {
    w <- stTraj[, c("target", "reference", "distance")]
  }

  # eps-neighborhood distance --------------------------------------------------

  # Assign NA to the trajectories for which all their states are out of the eps-neighborhood
  if (!is.null(eps)) {
    for (itarget in 1:length(targets_w)) {
      w[which(w$target == targets_w[itarget]), ]$distance <- ifelse(w[which(w$target == targets_w[itarget]), ]$reference %in% eps_neighborhood[[itarget]],
                                                                    w[which(w$target == targets_w[itarget]), ]$distance, NA)
    }
  }

  # k-neighbors ----------------------------------------------------------------

  # Assign NA to the trajectories for which all their states are not included in the k-nearest states
  for (itarget in 1:length(targets_w)) {
    w[which(w$target == targets_w[itarget]), ]$distance <- ifelse(w[which(w$target == targets_w[itarget]), ]$reference %in% k_neighbors[[itarget]],
                                                                  w[which(w$target == targets_w[itarget]), ]$distance, NA)
  }


  # Calculate weights ----------------------------------------------------------

  w <- data.table::data.table(w)

  # d_max
  d_max <- vapply(seq_along(targets_w), function(itarget){
    if (!is.null(eps) && !is.na(eps_w[itarget])) {
      max(eps_w[itarget], max(w[target == targets_w[itarget], ]$distance, na.rm = T))
    } else {
      max(w[target == targets_w[itarget], ]$distance, na.rm = T)
    }
  }, numeric(1))

  for (itarget in seq_along(targets_w)) {
    # NULL
    if (is.na(w_function_w[itarget])) {
      w[target == targets_w[itarget], w := 1]
    } else {
      # Linear
      if (w_function_w[itarget] == "linear") {
        w[target == targets_w[itarget], w := 1 - distance / d_max[itarget]]
      }

      # Power
      if (w_function_w[itarget] == "power") {
        w[target == targets_w[itarget], w := 1 - (distance / d_max[itarget])^alpha_w[itarget]]
      }

      # Exponential
      if (w_function_w[itarget] == "exponential") {
        w[target == targets_w[itarget], w := exp(-alpha_w[itarget] * distance / d_max[itarget])]
      }

      # Gaussian
      if (w_function_w[itarget] == "Gaussian") {
        w[target == targets_w[itarget], w := exp(-(distance / d_max[itarget])^2)]
      }

      # Hyperbolic
      if (w_function_w[itarget] == "hyperbolic") {
        w[target == targets_w[itarget], w := ((1 + distance / d_max[itarget])^(-alpha_w[itarget]) - 2^(-alpha_w[itarget])) / (1 - 2^(-alpha_w[itarget]))]
      }

      # Spherical
      if (w_function_w[itarget] == "spherical") {
        w[target == targets_w[itarget], w := 1 - 1.5*(distance / d_max[itarget]) + 0.5*(distance / d_max[itarget])^3]
      }
    }
  }

  w$w[is.na(w$w)] <- 0
  w <- data.frame(w)

  # Output ---------------------------------------------------------------------

  # Add data for the targets
  combined <- expand.grid(target = targets, reference = targets, w = 0)
  combined <- combined[which(!paste0(combined$target, "_", combined$reference) %in%
                               paste0(w$target, "_", w$reference)), ]

  # Targets with w_function = NA
  if (length(id_w) < length(targets)) {
    targets_w1 <- targets[-id_w]
    combined_w1 <- expand.grid(target = targets_w1,
                               reference = unique(trajectories[!trajectories %in% targets]),
                               w = 1)
    combined <- rbind(combined, combined_w1)
  }

  # Add targets data
  w <- rbind(w[, c("target", "reference", "w")], combined)

  # Repeat data for each trajectory state if it is going to be used in petra_edr
  if (w_per_state == T) {
    w <- merge(w, data.frame(reference = trajectories, states = states),
               by = "reference", all.x = T)
  }

  # Set the same order than in targets
  w <- lapply(setNames(targets, targets), function(itarget){
    dt <- w[w$target == itarget, ]
    dt[order(match(dt$reference, unique(trajectories))), ]$w
  })

  # Return a vector if there is only a target
  if (length(w) == 1) {
    w <- unlist(w, use.names = F)
  }

  return(w)

}

