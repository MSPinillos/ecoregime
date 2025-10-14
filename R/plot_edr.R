#' Plot Ecological Dynamic Regimes
#'
#' @description
#' Represents EDR trajectories in the state space. Trajectories and/or states can
#' be displayed in different colors based in a predefined classification or variable.
#'
#' @param x Symmetric matrix or `dist` object containing the dissimilarities
#' between each pair of states of all trajectories in the EDR. Alternatively, data
#' frame containing the coordinates of all trajectory states in an ordination space.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `x` belongs.
#' @param states Vector of integers indicating the order of the states in `x` for
#' each trajectory.
#' @param traj.colors Specification for the color of all individual trajectories
#' (defaults "grey") or a vector with length equal to the number of different
#' trajectories indicating the color for each individual trajectory.
#' @param state.colors Specification for the color of all trajectory states
#' (defaults equal to `traj.colors`), vector with length equal to the number
#' of states indicating the color for each trajectory state, or vector of colors
#' used to generate a gradient depending on the values of `variable` (if
#' `type = "gradient"`).
#' @param variable Numeric vector with equal length to the number of states to
#' be represented using a gradient of state colors (if `type = "gradient"`).
#' @param type One of the following `"trajectories"`, `"states"`, or `"gradient"`.
#' @param axes An integer vector indicating the pair of axes in the ordination
#' space to be plotted.
#' @param initial Flag indicating if the initial state must be plotted (only if
#' `type = "states"` or `type = "gradient"`)
#' @param ... Arguments for generic [plot()].
#'
#' @return
#' `plot_edr()` permits representing the trajectories of an Ecological Dynamic
#' Regime using different colors for each trajectory or state.
#'
#' @author Martina Sánchez-Pinillos
#'
#' @seealso
#' [`plot.RETRA()`] for plotting representative trajectories in an ordination space
#' representing the state space of the EDR.
#'
#' @export
#'
#' @examples
#'
#' # Data
#' state_variables <- EDR_data$EDR1$abundance
#' d <- EDR_data$EDR1$state_dissim
#'
#' # Coordinates in classic multidimensional scaling
#' x <- cmdscale(d, k = 3)
#'
#' # Plot trajectories 1-10 in "coral", 11-20 in "blue" and 21-30 in "gold"
#' plot_edr(x = x, trajectories = state_variables$traj,
#'          states = as.integer(state_variables$state),
#'          traj.colors = c(rep("coral", 10), rep("royalblue", 10), rep("gold", 10)),
#'          main = "type = 'trajectories'")
#' legend("bottomleft", legend = paste0("Trajectories ", c("1-10", "11-20", "21-30")),
#'        lty = 1, col = c("coral", "royalblue", "gold"))
#'
#' # Plot states with different colors depending on the state value
#' plot_edr(x = x, trajectories = state_variables$traj,
#'          states = as.integer(state_variables$state),
#'          traj.colors = NULL,
#'          state.colors = rep(RColorBrewer::brewer.pal(5, "Blues"),
#'                             length(unique(state_variables$traj))),
#'          type = "states", main = "type = 'states'")
#' legend("bottomleft", legend = paste0("State ", 1:5),
#'        pch = 15, col = RColorBrewer::brewer.pal(5, "Blues"))
#'
#' # Plot states with different colors depending on the abundance of sp1
#' plot_edr(x = x, trajectories = state_variables$traj,
#'          states = as.integer(state_variables$state),
#'          traj.colors = NULL, state.colors = viridis::viridis(5),
#'          variable = state_variables$sp1,
#'          type = "gradient", main = "type = 'gradient'", initial = TRUE)
#' legend("bottomleft",
#'        legend = c(paste0("abun sp1 = ", min(state_variables$sp1)),
#'                   rep(NA, 28),
#'                   paste0("abun sp1 = ", max(state_variables$sp1))),
#'        fill = viridis::viridis(30), border = NA, y.intersp = 0.2)


plot_edr <- function(x, trajectories, states, traj.colors = NULL, state.colors = NULL,
                     variable = NULL, type = "trajectories", axes = c(1, 2), initial = F, ...) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for x, trajectories, states
  if (!any(inherits(x, "matrix"),
           inherits(x, "data.frame"),
           inherits(x, "dist"))) {
    stop("'x' needs to be of any of these classes: 'matrix', 'data.frame', 'dist'.")
  }
  if (any(inherits(x, "matrix"), inherits(x, "data.frame")) &
      nrow(x) != length(trajectories)) {
    stop("The length of 'trajectories' must be equal to the number of rows of 'x'.")
  } else if (inherits(x, "dist") & nrow(as.matrix(x)) != length(trajectories)) {
    stop("The length of 'trajectories' must be equal to the number of rows of 'x'.")
  }
  if (length(trajectories) != length(states)) {
    stop("The length of 'trajectories' must be equal to the length of 'states'.")
  }
  if (!is.integer(states)) {
    stop("'states' needs to be of class integer.")
  }

  # Check type
  type <- match.arg(type, c("trajectories", "states", "gradient"))

  # Check and scale variable
  if (type == "gradient") {
    if (any(!is.numeric(variable), is.null(variable), length(variable) != length(states))) {
      stop(cat("'variable' must be a numeric vector of length", length(states)))
    }
    if (any(is.na(variable))) {
      stop("There are missing values in 'variable'")
    }
    sc_var <- round(100 * (variable - min(variable)) / (max(variable) - min(variable)) + 1)
  }

  # Set order
  traj_st <- data.frame(trajectories = trajectories, states = states)
  set_order <- order(traj_st$trajectories, traj_st$states)

  # STATE COORDINATES ----------------------------------------------------------

  # Coordinates MDS
  if (inherits(x, "dist") || isSymmetric(as.matrix(x))) {
    warning(cat("Trajectories will be displayed in an ordination space generated through multidimensional scaling (MDS). You can avoid this step by providing state coordinates in 'x'.", "\n"))
    st_coord <- data.frame(smacof::mds(delta = x, ndim = ncol(x)-1,
                                       itmax = 300, verbose = F)$conf)
  } else {
    st_coord <- data.frame(x)
  }

  if (any(set_order != 1:length(trajectories))) {
    st_coord <- st_coord[set_order, ]
    st_coord.ls <- lapply(unique(trajectories), function(itraj){
      dt <- st_coord[which(sort(trajectories) == itraj), ]
    })
  } else {
    st_coord.ls <- lapply(unique(trajectories), function(itraj){
      dt <- st_coord[which(trajectories == itraj), ]
    })
  }

  # DEFINE COLORS --------------------------------------------------------------

  # Define trajectory colors
  if (!is.null(traj.colors)) {
    if (length(traj.colors) == 1) {
      traj.colors = rep(traj.colors, length(unique(trajectories)))
    }
    if (length(traj.colors) != length(unique(trajectories))) {
      warning("The length of 'traj.colors' is different from the number of trajectories. Only the first element will be used.")
      traj.colors = rep(traj.colors[1], length(unique(trajectories)))
    }
  }
  if (is.null(traj.colors)) {
    traj.colors = rep("grey", length(unique(trajectories)))
  }

  # Define state colors
  if (!is.null(state.colors)) {
    if (length(state.colors) == 1) {
      state.colors = rep(state.colors, length(states))
    }
    if (type == "states" & length(state.colors) < length(states)) {
      warning("'state.colors' has a shorter length of 'states'. Only the first element will be used.")
      state.colors = rep(state.colors[1], length(states))
    }
    if (type == "gradient") {
      colorfunc <- grDevices::colorRampPalette(state.colors)
      state.colors <- colorfunc(101)[sc_var]
    }
  }
  if (any(set_order != 1:length(trajectories))) {
    state.colors <- state.colors[set_order]
  }

  if (any(set_order != 1:length(trajectories))) {
    state.colors.ls <- lapply(unique(trajectories), function(itraj){
      state.colors[which(sort(trajectories) == itraj)]
    })
  } else {
    state.colors.ls <- lapply(unique(trajectories), function(itraj){
      state.colors[which(trajectories == itraj)]
    })
  }

  if (is.null(state.colors)) {
    state.colors.ls = lapply(seq_along(unique(trajectories)), function(itraj){
      rep(traj.colors[itraj], sum(trajectories == unique(trajectories)[itraj]))
    })
  }

  # PLOT EDR ----------------------------------------------------------

  # Plot individual trajectories in the EDR
  plot(st_coord[, axes], type = "n",
       xlab = paste0("Axis ", axes[1]),
       ylab = paste0("Axis ", axes[2]), ...)

  if (type == "trajectories") {
    for (itraj in seq_along(st_coord.ls)) {
      istate = 1
      while (istate < nrow(st_coord.ls[[itraj]])) {
        shape::Arrows(x0 = st_coord.ls[[itraj]][istate, axes[1]], y0 = st_coord.ls[[itraj]][istate, axes[2]],
                      x1 = st_coord.ls[[itraj]][istate + 1, axes[1]], y1 = st_coord.ls[[itraj]][istate + 1, axes[2]],
                      col = traj.colors[itraj], arr.adj = 1)
        istate = istate + 1
      }
    }
  } else if (type %in% c("states", "gradient")) {
    for (itraj in seq_along(st_coord.ls)) {
      istate = 1
      if (any(initial == T, nrow(st_coord.ls[[itraj]]) == 1)) {
        graphics::points(x = st_coord.ls[[itraj]][istate, axes[1]], y = st_coord.ls[[itraj]][istate, axes[2]],
                         col = state.colors.ls[[itraj]][istate], pch = 20)
      }
      while (istate < nrow(st_coord.ls[[itraj]])) {
        shape::Arrows(x0 = st_coord.ls[[itraj]][istate, axes[1]], y0 = st_coord.ls[[itraj]][istate, axes[2]],
                      x1 = st_coord.ls[[itraj]][istate + 1, axes[1]], y1 = st_coord.ls[[itraj]][istate + 1, axes[2]],
                      col = traj.colors[itraj], arr.adj = 1, arr.type = "circle",
                      arr.col = state.colors.ls[[itraj]][istate], arr.length = 0.01)

        shape::Arrows(x0 = st_coord.ls[[itraj]][istate, axes[1]], y0 = st_coord.ls[[itraj]][istate, axes[2]],
                      x1 = st_coord.ls[[itraj]][istate + 1, axes[1]], y1 = st_coord.ls[[itraj]][istate + 1, axes[2]],
                      col = traj.colors[itraj], arr.adj = 1, lty = 0,
                      arr.col = state.colors.ls[[itraj]][istate+1])
        istate = istate + 1
      }
    }
  }
}
