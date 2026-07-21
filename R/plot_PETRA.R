#' Plot predicted trajectories in Ecological Dynamic Regimes
#'
#' @description
#' Represent predicted trajectories in the state space of an Ecological Dynamic Regime
#' (EDR) distinguishing between the observed states belonging to the target
#' and the forecasted states using the function `petra_edr()`.
#'
#' @param x Object of class `PETRA` returned by `petra_edr()` using `return_args = T`.
#' @param petra.colors Specification for the color of all predicted trajectories
#' (defaults red) or a vector with length equal to the number of predicted
#' trajectories in `x` indicating the color for each predicted trajectory.
#' @param target.colors Specification for the color of the observed states (target)
#' within a predicted trajectory in `x` or a vector with length equal to the number
#' of predicted trajectories in `x` indicating a color for each target.
#' @param traj.colors Specification for the color of the observed trajectories in
#' the EDR (defaults grey).
#' @param uncert.metric Column name of `x$predicted_dist` used to display predicted
#' states with a gradient color depending on that metric. One of `"N"`, `"mean_dist"`,
#' `"min_dist"`, `"max_dist"`.
#' @param uncert.colors Vector of colors used to generate a gradient depending on
#' the values of `uncert.metric`.
#' @param uncert.range Numeric vector including the minimum and maximum values to
#' scale `uncert.metric`.
#' @param coord Data frame containing the coordinates of all trajectory states in
#' an ordination space (including the predicted states).
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `coord` belongs.
#' @param states Vector of integers indicating the order of the states in `coord`
#' for each trajectory.
#' @param axes An integer vector indicating the pair of axes in the ordination
#' space to be plotted.
#' @param ... Arguments for [`plot_edr()`].
#'
#' @return
#' Graphical representation of a set of individual trajectories and the predicted
#' trajectories in an ordination space defined by `coord` or calculated by applying
#' metric multidimensional scaling (mMDS; Borg and Groenen, 2005) to a dissimilarity
#' matrix computed from the arguments used in `petra_edr()` to generate `x`.
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.).
#' Springer.
#'
#' Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
#' Forecasting ecological trajectories from ecological dynamic regimes to improve
#' resilience analysis. *Methods in Ecology and Evolution*. <doi:10.1111/2041-210x.70372>
#'
#' @seealso
#' [`petra_edr()`] for predicting ecological trajectories from EDRs.
#'
#' [`MPD()`] for estimating the prediction accuracy of `petra-edr()` outputs.
#'
#' [`plot_edr()`] for representing EDRs in a multidimensional state space.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'
#' # Compute PETRA-EDR ---------------------------------------------------------
#'
#' # State variables of the states in the EDR and the target
#' EDR_var <- EDR_data$EDR3$abundance
#' target_var <- EDR_data$EDR3_disturbed$abundance[disturbed_states == 0]
#'
#' # Define state_var including the state variables for the states in the EDR and
#' # the targets
#' state_var <- data.frame(rbind(EDR_var[, paste0("sp", 1:12)],
#'                               target_var[, paste0("sp", 1:12)]))
#'
#' # Compute PETRA-EDR for each target. Set return_args = TRUE
#' petra <- petra_edr(state_var = state_var,
#'                    trajectories = c(EDR_var$traj, target_var$traj),
#'                    states = c(EDR_var$state, target_var$state),
#'                    targets = unique(target_var$traj),
#'                    k = 5L, minPts = 2L,
#'                    d_function = "vegan::vegdist",
#'                    d_args = list(x = state_var, method = "bray"),
#'                    return_args = TRUE)
#'
#' # Recalculate state_var and d, including the data for the predicted trajectories
#' state_var <- rbind(EDR_var[, paste0("sp", 1:12)],
#'                    petra$state_var)
#'
#' d <- vegan::vegdist(state_var, method = "bray")
#'
#' # Compute PCoA (optional)
#' pcoa <- cmdscale(d = d)#'
#'
#' # Example 1 -----------------------------------------------------------------
#'
#' # Display each predicted trajectories with different colors
#' plot(x = petra,
#'      petra.colors = c("red", "royalblue", "seagreen"),
#'      traj.colors = "grey",
#'      coord = pcoa,
#'      trajectories = c(EDR_var$traj, petra$trajectories),
#'      states = as.integer(c(EDR_var$state, petra$states)),
#'      xlab = "MDS D1", ylab = "MDS D2",
#'      main = "Predicted trajectories in EDR")
#' legend("bottomleft", legend = unique(target_var$traj),
#'        lwd = 2,
#'        col = c("red", "royalblue", "seagreen"),
#'        title = "Target")
#'
#' # Example 2 -----------------------------------------------------------------
#'
#' # Identify the observed states with a different color
#' plot(x = petra,
#'      petra.colors = "black",
#'      target.colors = c("red", "royalblue", "seagreen"),
#'      traj.colors = "grey",
#'      coord = pcoa,
#'      trajectories = c(EDR_var$traj, petra$trajectories),
#'      states = as.integer(c(EDR_var$state, petra$states)),
#'      xlab = "MDS D1", ylab = "MDS D2",
#'      main = "Observed and predicted states")
#' legend("bottomleft",
#'        legend = c(paste0("Target ", unique(target_var$traj)), "Predicted states"),
#'        lwd = 2,
#'        col = c("red", "royalblue", "seagreen", "black"))
#'
#' # Example 3 -----------------------------------------------------------------
#'
#' # Display predicted states based on an uncertainty metric
#' plot(x = petra,
#'      petra.colors = "black",
#'      target.colors = "black",
#'      traj.colors = "grey",
#'      uncert.metric = "N",
#'      uncert.colors = hcl.colors(4, rev = TRUE),
#'      coord = pcoa,
#'      trajectories = c(EDR_var$traj, petra$trajectories),
#'      states = as.integer(c(EDR_var$state, petra$states)),
#'      xlab = "MDS D1", ylab = "MDS D2",
#'      main = "Uncertainty of the predicted states")
#' legend("bottomleft", legend = c(paste0("N = ", min(petra$predicted_dist$N)),
#'                                 rep(NA, 18),
#'                                 paste0("N = ", max(petra$predicted_dist$N))),
#'        fill = hcl.colors(20, rev = TRUE), border = NA, y.intersp = 0.2)
#' }
#'
#'

plot.PETRA <- function(x, petra.colors = "red", target.colors = NULL, traj.colors = "grey",
                       uncert.metric = NULL, uncert.colors = NULL, uncert.range = NULL,
                       coord = NULL, trajectories = NULL, states = NULL, axes = c(1, 2), ...) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format of x
  if (!"arguments" %in% names(x)) {
    stop("'x' must contain 'arguments'. Set `return_arguments = T` when computing 'petra_edr()' to get 'x'.")
  }

  # Check uncert.metric, uncert.colors, uncert.range
  if (!is.null(uncert.metric)) {
    # Check uncert.metric
    uncert.metric <- match.arg(uncert.metric, names(x$predicted_dist)[-c(1:2)])
    # Check uncert.colors
    if (length(uncert.colors) < 2) {
      stop("Provide at least two non-NA values in 'uncert.colors'.")
    }
    # Check uncert.range
    if (!is.null(uncert.range)) {
      if (!is.numeric(uncert.range)) {
        stop("Invalid 'uncert.range' value.")
      }
      if (min(uncert.range) > min(x$predicted_dist[[uncert.metric]])) {
        stop(paste0("The minimum value in 'uncert.range' cannot be greater than ", min(x$predicted_dist[[uncert.metric]])))
      }
      if (max(uncert.range) < max(x$predicted_dist[[uncert.metric]])) {
        stop(paste0("The maximum value in 'uncert.range' cannot be lower than ", max(x$predicted_dist[[uncert.metric]])))
      }
    }
  }

  # Check coordinates
  if (!is.null(coord)) {
    if (!any(inherits(coord, "matrix"),
             inherits(coord, "data.frame"))) {
      stop("'coord' needs to be of any of these classes: 'matrix', 'data.frame'.")
    }
    if (any(is.null(trajectories), is.null(states))) {
      stop("Provide values for 'trajectories' and 'states'.")
    }
    if (any(nrow(coord) != length(trajectories), nrow(coord) != length(states))) {
      stop("The length of 'trajectories' and 'states' must be equal to the number of rows of 'coord'")
    }
    if (!is.integer(states)) {
      stop("'states' needs to be of class integer.")
    }
    if (!all(c(paste0(x$arguments$trajectories, "_", x$arguments$states),
               paste0(x$trajectories, "_", x$states)) %in% paste0(trajectories, "_", states))) {
      stop("'coord' must contain coordinates for all states in the EDR and the predicted trajectories in 'x'")
    }
  }

  ## DEFINE VARIABLES ----------------------------------------------------------

  # Variables used to compute x
  arg_state_var <- x$arguments$d_args[[x$arguments$args_state_var]]
  arg_trajectories <- x$arguments$trajectories
  arg_states <- x$arguments$states
  arg_targets <- x$arguments$targets

  # ID of the targets in the original data
  ID_target <- which(arg_trajectories %in% arg_targets)

  # Remove the targets from the original data and add predicted trajectories
  new_trajectories <- c(arg_trajectories[-ID_target], x$trajectories)
  new_states <- as.integer(c(arg_states[-ID_target], x$states))

  # Check if the states in coord is in order and PETRAs are at the end
  if (!is.null(coord)) {
    new_traj_st <- paste0(new_trajectories, "_", new_states)
    traj_st <- paste0(trajectories, "_", states)

    # Set the right order if needed (new_trajectories and new_states)
    if (identical(new_traj_st, traj_st) == F) {
      set_order <- match(new_traj_st, traj_st)
      if (length(traj_st) > length(new_traj_st)) {
        warning("Only the trajectories used to compute 'x' will be displayed. Use plot_edr() to plot additional trajectories.")
      }
      coord <- coord[set_order[!is.na(set_order)], ]
    }
  }

  # Compute a dissimilarity matrix and MDS
  if (is.null(coord)) {
    warning("Trajectories will be displayed in an ordination space generated through multidimensional scaling (MDS). You can avoid this step by providing state coordinates in 'coord'.")

    if (inherits(arg_state_var, "list")) {
      new_state_var <- c(arg_state_var[-ID_target], x$state_var)
      class(new_state_var) <- class(arg_state_var)
    } else {
      new_state_var <- rbind(arg_state_var[-ID_target, ], x$state_var)
    }

    # Identify the argument in d_function corresponding to state_var
    istate_var <- x$arguments$args_state_var

    # Dissimilarity between the states of the predicted trajectories and each state of the EDR
    new_d <- as.matrix(x$arguments$d)[-ID_target, -ID_target]
    d_function <- x$arguments$d_function
    d_args <- x$arguments$d_args
    d_pre_obs <- numeric()

    for(ipredicted in which(new_trajectories %in% x$trajectories)) {
      # Calculate the dissimilarity of each predicted state and the observed states in the EDR
      for (iobserved in 1:ipredicted) {
        if (inherits(new_state_var, "list")) {
          d_args[[istate_var]] <- c(new_state_var[iobserved], new_state_var[ipredicted])
          class(d_args[[istate_var]]) <- class(arg_state_var)
          d_pre_obs <- c(d_pre_obs, do.call(eval(parse(text = d_function)), args = d_args))
        } else {
          d_args[[istate_var]] <- rbind(new_state_var[iobserved, ], new_state_var[ipredicted, ])
          d_pre_obs <- c(d_pre_obs, do.call(eval(parse(text = d_function)), args = d_args))
        }
      }
      # Add the values to the dissimilarity matrix
      new_d <- rbind(new_d, d_pre_obs[1:ncol(new_d)])
      new_d <- cbind(new_d, d_pre_obs)
      row.names(new_d)[ipredicted] <- paste0(new_trajectories[ipredicted], new_states[ipredicted])
      colnames(new_d)[ipredicted] <- paste0(new_trajectories[ipredicted], new_states[ipredicted])

      # Reset d_pre_obs for the next loop
      d_pre_obs <- numeric()
    }

    # Apply MDS
    coord <- data.frame(smacof::mds(delta = new_d, ndim = ncol(new_d)-1,
                                    itmax = 300, verbose = F)$conf)
  }


  ## COLORS --------------------------------------------------------------------

  # Define a color for each trajectory in the EDR
  if (length(traj.colors) == 1) {
    traj.colors <- rep(traj.colors, length(unique(arg_trajectories[-ID_target])))
  } else if (length(traj.colors) != length(unique(arg_trajectories[-ID_target]))) {
    warning("The length of 'traj.colors' must be 1 or equal to the number of trajectories in the EDR. \nOnly the first element of 'traj.colors' will be used.")
    traj.colors <- rep(traj.colors[1], length(unique(arg_trajectories[-ID_target])))
  }

  # Define a color for each predicted trajectory
  if (length(petra.colors) == 1) {
    petra.colors <- rep(petra.colors, length(unique(x$trajectories)))
  }
  if (length(petra.colors) != length(unique(x$trajectories))) {
    warning("The length of 'petra.colors' is different from the number of predicted trajectories. Only the first element of 'petra.colors' will be used.")
    petra.colors <- rep(petra.colors[1], length(unique(x$trajectories)))
  }

  if (any(!is.null(target.colors), !is.null(uncert.metric))) {
    # Trajectory-state in the arguments of PETRA and the new trajectories
    # including the predicted states
    arg_traj_st <- paste0(arg_trajectories, "_", arg_states)
    new_traj_st <- paste0(new_trajectories, "_", new_states)

    # Color of the EDR and PETRA states

    # edr_states <- rep(unique(traj.colors), length(arg_trajectories[-ID_target]))
    edr_states <- data.table::merge.data.table(data.table::data.table(traj = unique(arg_trajectories[-ID_target]),
                                                                      traj.colors = traj.colors),
                                               data.table::data.table(traj = arg_trajectories[-ID_target]),
                                               sort = FALSE)$traj.colors


    petra_states <- unlist(lapply(seq_along(unique(x$trajectories)), function(itarget){
      rep(petra.colors[itarget],
          sum(x$trajectories == unique(x$trajectories)[itarget]))
    }))
    state.colors <- c(edr_states, petra_states)
  }


  # Define the target states and their colors
  if (!is.null(target.colors)) {

    # Identify observed states
    observed_states <- which(new_traj_st %in% arg_traj_st[ID_target])

    # Index of the observed states
    petra_observed <- lapply(seq_along(unique(x$trajectories)), function(ipetra){
      ipetra_traj <- unique(x$trajectories)[ipetra]
      ipetra_states <- which(new_traj_st %in% paste0(x$trajectories, "_", x$states)[x$trajectories == ipetra_traj])
      ipetra_observed <- observed_states[which(observed_states %in% ipetra_states)]
      return(ipetra_observed)
    })

    # Color of the target states
    if (length(target.colors) == 1) {
      state.colors[observed_states] <- target.colors
      target.colors <- rep(target.colors, length(unique(x$trajectories)))
    }
    if (length(target.colors) == length(unique(x$trajectories))) {
      for (ipetra in seq_along(unique(x$trajectories))) {
        state.colors[petra_observed[[ipetra]]] <- target.colors[ipetra]
      }
    }
  }

  if (!is.null(uncert.metric)) {

    # Index of the predicted states
    predicted_states <- which(!new_traj_st %in% arg_traj_st)

    # Scale the values of uncert.metric
    if (!is.null(uncert.range)) {
      uncert.values <- c(x$predicted_dist[[uncert.metric]], uncert.range)
    } else {
      uncert.values <- x$predicted_dist[[uncert.metric]]
    }

    sc_var <- round(100 * (uncert.values - min(uncert.values)) / (max(uncert.values) - min(uncert.values)) + 1)

    # Define colors of the predicted states
    colorfunc <- colorRampPalette(uncert.colors)
    predicted.colors <- colorfunc(101)[sc_var]

    # Replace colors for predicted states
    state.colors[predicted_states] <- predicted.colors[1:length(x$predicted_dist[[uncert.metric]])]

    # Index of the first petra state in new_trajectories
    first_petra_st <- vapply(unique(x$trajectories), function(ipetra){
      ipetra_traj <- which(x$trajectories == ipetra)
      which(new_traj_st %in% paste0(x$trajectories[ipetra_traj], "_", x$states[ipetra_traj])[1])
    }, integer(1))

  }


  ## APPLY plot_edr() ----------------------------------------------------------

  if (is.null(target.colors)) {
    plot_edr(x = coord, trajectories = new_trajectories, states = new_states,
             traj.colors = c(traj.colors, petra.colors), ...)
  }

  if (any(!is.null(target.colors), !is.null(uncert.metric))) {
    plot_edr(x = coord, trajectories = new_trajectories, states = new_states,
             traj.colors = c(traj.colors, petra.colors),
             state.colors = state.colors, type = "states", ...)
    if (!is.null(target.colors)) {
      for (ipetra in seq_along(unique(x$trajectories))) {
        lines(x = coord[petra_observed[[ipetra]], axes[1]],
              y = coord[petra_observed[[ipetra]], axes[2]], col = target.colors[ipetra])
        if (which(new_trajectories == unique(x$trajectories)[[ipetra]])[1] == petra_observed[[ipetra]][1]) {
          points(x = coord[petra_observed[[ipetra]][1], axes[1]],
                 y = coord[petra_observed[[ipetra]][1], axes[2]],
                 col = state.colors[petra_observed[[ipetra]][1]], pch = 20)
        }
      }
    }
    if (!is.null(uncert.metric)) {
      points(x = coord[first_petra_st, axes[1]], y = coord[first_petra_st, axes[2]],
             col = state.colors[first_petra_st], pch = 20)
    }
  }
}

