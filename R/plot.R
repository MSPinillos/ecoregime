#' Plot representative trajectories of Ecological Dynamic Regimes
#'
#' @description
#' Plot representative trajectories of an Ecological Dynamic Regime (EDR) in the
#' state space, distinguishing between the segments belonging to real trajectories
#' of the EDR and the artificial links between segments.
#'
#' @param x Object of class `RETRA`.
#' @param d Symmetric matrix or `dist` object containing the dissimilarities
#' between each pair of states of all trajectories in the EDR or data frame
#' containing the coordinates of all trajectory states in an ordination space.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `d` belongs.
#' @param states Vector of integers indicating the order of the states in `d` for
#' each trajectory.
#' @param select_RT Optional string indicating the name of a representative
#' trajectory that must be highlighted in the plot. By default (`select_RT` = `NULL`),
#' all representative trajectories are represented with the same color.
#' @param traj.colors Specification for the color of all individual trajectories
#' (defaults "grey") or a vector with length equal to the number of trajectories
#' indicating the color for each individual trajectory.
#' @param RT.colors Specification for the color of representative trajectories
#' (defaults "black").
#' @param sel.color Specification for the color of the selected representative
#' trajectory (defaults "red"). Only if `!is.null(select_RT)`.
#' @param link.color Specification for the color of the links between trajectory
#' segments forming representative trajectories. By default, the same color than
#' `RT.colors` is used.
#' @param link.lty The line type of the links between trajectory segments forming
#' representative trajectories. Defaults 2 = dashed (See [graphics::par]).
#' @param axes An integer vector indicating the pair of axes in the ordination
#' space to be plotted.
#' @param ... Arguments for generic [plot()].
#'
#' @return
#' The function `plot()` plots a set of individual trajectories and the
#' representative trajectories in an ordination space defined through `d` or
#' calculated by applying metric multidimensional scaling (mMDS; Borg and Groenen,
#' 2005) to `d`.
#'
#' @author Martina Sánchez-Pinillos, CNRS, Univ. Montpellier
#'
#' @references
#' Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. Ecological Dynamic
#' Regimes: Identification, characterization, and comparison
#'
#' De Cáceres, M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ,
#' Condit R & Hubbell S. (2019). Trajectory analysis in community ecology. Ecological
#' Monographs.
#'
#' Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling (2nd ed.).
#' Springer.
#'
#' @note
#' This function uses a modified version of the [ecotraj::trajectoryPlot()] function
#' in 'ecotraj' (v.0.0.3; De Cáceres et al. 2019). The modification was done the
#' 2022-12-18. It allows using the function [shape::Arrows()] to plot ecological
#' trajectories instead of [graphics::arrows()] and setting the parameters of
#' generic [plot()].
#'
#' @seealso
#' [`retra_edr()`] for identifying representative trajectories in EDRs applying
#' RETRA-EDR.
#'
#' [`define_retra()`] for defining representative trajectories from a subset of
#' segments or trajectory features.
#'
#' [`summary()`] for summarizing representative trajectories in EDRs.
#'
#' @export
#'
#' @examples
#' # Example 1 -----------------------------------------------------------------
#'
#' # d contains the dissimilarities between trajectory states
#' d <- EDR_data$EDR1$state_dissim
#'
#' # trajectories and states are defined according to `d` entries.
#' trajectories <- EDR_data$EDR1$abundance$traj
#' states <- EDR_data$EDR1$abundance$state
#'
#' # x defined from retra_edr(). We obtain three representative trajectories.
#' RT <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)
#' summary(RT)
#'
#' # Plot individual trajectories in blue and representative trajectories in orange,
#' # "T2" will be displayed in green. Artificial links will be displayed with a
#' # dotted line.
#' plot(x = RT, d = d, trajectories = trajectories, states = states, select_RT = "T2",
#'      traj.colors = "lightblue", RT.colors = "orange", sel.color = "darkgreen",
#'      link.lty = 3, main = "Representative trajectories in EDR1")
#'
#' # Example 2 -----------------------------------------------------------------
#'
#' # d contains the coordinates in an ordination space. For example, we use
#' # the coordinates of the trajectory states after applying a principal component
#' # analysis (PCA) to an abundance matrix.
#' abun <- EDR_data$EDR1$abundance
#' pca <- prcomp(abun[, -c(1:3)])
#' coord <- data.frame(pca$x)
#'
#' # trajectories and states are defined according to the abundance matrix
#' # used in the PCA
#' trajectories <- EDR_data$EDR1$abundance$traj
#' states <- EDR_data$EDR1$abundance$state
#'
#' # Instead of using the representative trajectories obtained from `retra_edr()`,
#' # we will define the set of trajectories that we want to highlight. For example,
#' # we can select the trajectories whose initial and final states are in the
#' # extremes of the first axis.
#' T1 <- trajectories[which.max(coord[, 1])]
#' T2 <- trajectories[which.min(coord[, 1])]
#' RT_traj <- c(trajectories[trajectories %in% T1],
#'              trajectories[trajectories %in% T2])
#' RT_states <- c(states[which(trajectories %in% T1)],
#'                states[which(trajectories %in% T2)])
#'
#' # Create a data frame to generate a RETRA object using define_retra
#' RT_df <- data.frame(RT = c(rep("T1", sum(trajectories %in% T1)),
#'                            rep("T2", sum(trajectories %in% T2))),
#'                  RT_traj = RT_traj,
#'                  RT_states = as.integer(RT_states))
#' RT_retra <- define_retra(data = RT_df)
#'
#' # Plot the defined trajectories with the default graphic values
#' plot(x = RT_retra, d = coord, trajectories = trajectories, states = states,
#'      main = "Extreme trajectories in EDR1")
#'

plot.RETRA <- function (x, d, trajectories, states, select_RT = NULL,
                        traj.colors = NULL, RT.colors = NULL, sel.color = NULL,
                        link.color = NULL, link.lty = 2, axes = c(1, 2), ...) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the formats
  if (!is.data.frame(d)) {
    if (all(!is.matrix(d), !inherits(d, "dist")) |
        nrow(as.matrix(d)) != ncol(as.matrix(d))) {
      stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist' containing state dissimilarities. Alternatively, you can use a data frame containing state coordinates in a ordination space.")
    }
  }
  if (is.data.frame(d) & sum(unlist(lapply(d, is.numeric)) == F) > 0) {
    stop("All entries in 'd' need to be numeric.")
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

  ## RETRESENTATIVE TRAJECTORIES -----------------------------------------------

  # Number of trajectories
  nRT <- length(x)

  RT_names <- character(0)
  RT_traj <- character(0)
  RT_states <- integer(0)
  for (iRT in seq_len(nRT)) {
    # Representative segments
    segs <- x[[iRT]]$Segments
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segs)), "-")

    # Representative trajectories based on the number of states
    RT_names <- c(RT_names, rep(names(x)[iRT], 2*length(segs)))
    for (iseg in seg_components) {
      RT_traj <- c(RT_traj, rep(iseg[1], 2))
      RT_states <- as.integer(c(RT_states, iseg[2:3]))
    }
  }

  if(sum(!paste0(RT_traj, RT_states) %in% paste0(trajectories, states)) > 0){
    stop("The states in forming representative trajectories (identified through 'Segments' in x) need to be included in 'trajectories' and 'states'.")
  }


  if (!is.null(select_RT) && !(select_RT %in% RT_names)) {
    stop(cat("select_RT = \"", select_RT, "\"; \"", select_RT, "\" needs to be included in 'x' or 'x$RT_names'", sep = ""))
  }


  # STATE COORDINATES ----------------------------------------------------------

    # Coordinates MDS
  if (inherits(d, "dist") || isSymmetric(as.matrix(d))) {
    warning(cat("Representative trajectories will be displayed in an ordination space generated through multidimensional scaling (MDS). You can avoid this step by providing state coordinates in the 'd' argument.", "\n"))
    statesMDS <- data.frame(smacof::mds(delta = d, ndim = ncol(d)-1,
                                        itmax = 300, verbose = F)$conf)
  } else {
    warning(cat("Representative trajectories will be displayed considering the coordinates provided in 'd'."))
    statesMDS <- data.frame(d)
  }
  statesMDS$ID <- paste0(trajectories, "_", states)

  # RT data: ID_plot, RT_ID
  RT_data <- data.frame(ID = paste0(RT_traj, "_", RT_states), RT_traj = RT_traj, RT_states = RT_states,
                        RT = RT_names, order_RT_states = unlist(lapply(unique(RT_names), function(iT){
                          seq_along(RT_names[RT_names == iT])
                        })))
  ID_RT <- names(sort(table(RT_names)))
  if (nRT > 1) {
    # Place the selected RT to be displayed at the end
    if (!is.null(select_RT)) {
      ID_RT <- c(ID_RT[-which(ID_RT == select_RT)], select_RT)
    }
    RT_data <- RT_data[order(match(RT_data$RT, ID_RT)), ]
    RT_traj <- RT_data$RT_traj
    RT_states <- RT_data$RT_states
  }

  # Coordinates of RT states
  RT_coords <- merge(RT_data, statesMDS, by = "ID", all.d = T)
  RT_coords <- RT_coords[match(paste0(RT_data$RT, "_", RT_data$order_RT_states),
                               paste0(RT_coords$RT, "_", RT_coords$order_RT_states)), ]

  links <- lapply(1:nRT, function(iT){
    ind_traj <- which(RT_data$RT == ID_RT[iT])

    nStates <- nrow(RT_data[ind_traj, ])
    ind_links <- integer()
    i <- 2
    while (i <= nStates) {
      if (!all(RT_data[ind_traj, ]$RT_traj[i] == RT_data[ind_traj, ]$RT_traj[(i-1)],
               RT_data[ind_traj, ]$RT_states[i] - RT_data[ind_traj, ]$RT_states[(i-1)] <= 1)) {
        ind_links <- c(ind_links, i-1, i)
      }
      i <- i+1
    }

    links <- RT_data[ind_traj, ][ind_links, ]

  })
  links <- data.frame(data.table::rbindlist(links))

  if (nrow(links) > 0) {
    # Add and ID for the link and the order of the states
    links$Link <- rep(1:(nrow(links)/2), each = 2)
    links$Link_state <- rep(1:2, nrow(links)/2)
    links <- links[, c("ID", "RT", "Link", "Link_state")]
    # Coordinates of link states
    Lk_coords <- merge(links, statesMDS, by = "ID", all.d = T)
    Lk_coords <- Lk_coords[match(paste0(links$RT, "_", links$Link, "_", links$Link_state),
                                 paste0(Lk_coords$RT, "_", Lk_coords$Link, "_", Lk_coords$Link_state)), ]
  }

  # PLOT TRAJECTORIES ----------------------------------------------------------

  # Plot all trajectories in the EDR
  if (!is.null(traj.colors)) {
    if (length(traj.colors) == 1) {
      traj.colors = rep(traj.colors, length(unique(trajectories)))
    }
    if (length(traj.colors) < length(unique(trajectories))) {
      warning("'traj.colors' has a shorter length than the number of trajectories and only the first element will be used.")
      traj.colors = rep(traj.colors[1], length(unique(trajectories)))
    }
  }
  if (is.null(traj.colors)) {
    traj.colors = rep("grey", length(unique(trajectories)))
  }
  trajectoryPlot2(d = statesMDS[, -ncol(statesMDS)], sites = trajectories,
                  surveys = states, traj.colors = traj.colors, axes = axes, ...)

  # Plot segments
  if (is.null(RT.colors)) {
    RT.colors <- rep("black", nRT)
  }
  if (length(RT.colors) == 1) {
    RT.colors <- rep(RT.colors, nRT)
  }
  if (!is.null(select_RT)) {
    if (!is.null(sel.color)) {
      RT.colors[nRT] <- sel.color
    }
    if (is.null(sel.color)) {
      RT.colors[nRT] <- "red"
    }
  }

  for(iT in seq_len(nRT)){
    ind_RT <- which(RT_data$RT == ID_RT[iT])
    ind_dup <- which(duplicated(RT_data[ind_RT, 1:4]))
    if (length(ind_dup) > 0) {
      ind_RT <- ind_RT[-ind_dup]
    }

    xp <- RT_coords[ind_RT, 5+axes[1]]
    yp <- RT_coords[ind_RT, 5+axes[2]]

    iT_traj <- RT_data[ind_RT, ]$RT_traj
    iT_states <- RT_data[ind_RT, ]$RT_states
    ID_traj <- unique(iT_traj)

    for (i in 1:length(ID_traj)) {
      ind_states = which(iT_traj == ID_traj[i])
      ind_states = ind_states[order(iT_states[iT_traj == ID_traj[i]])]
      if (length(ind_states) > 1) {
        for (t in 1:(length(ind_states) - 1)) {
          niini = ind_states[t]
          nifin = ind_states[t + 1]
          if (nifin - niini == 1) {
            shape::Arrows(xp[niini], yp[niini], xp[nifin], yp[nifin],
                          col = RT.colors[iT], arr.adj = 1)
          }
        }
      }
    }
  }

  # Plot the links
  if (!is.null(link.color)) {
    link.color <- rep(link.color, nRT)
  } else {
    link.color <- RT.colors
  }
  if (is.null(link.lty)){
    link.lty <- 2
  }
  for(iT in seq_len(nRT)){
    ind_RT <- which(links$RT == ID_RT[iT])

    if (length(ind_RT) > 0) {
      xp <- Lk_coords[ind_RT, 4+axes[1]]
      yp <- Lk_coords[ind_RT, 4+axes[2]]

      iL_traj <- links[ind_RT, ]$Link
      ID_link <- unique(iL_traj)

      for (i in 1:length(ID_link)) {
        ind_states = which(iL_traj == ID_link[i])
        lines(x = xp[ind_states], y = yp[ind_states], col = link.color[iT], lty = link.lty)
      }
    }
  }
}



# Function to plot individual trajectories. This is a modified version of the
# function trajectoryPlot() in 'ecotraj' Version 0.0.3 (De Cáceres et al. 2019).
# The modification allows plotting ecological trajectories as arrows using the
# Arrows() function in 'shape', instead of the arrows() function in base.
trajectoryPlot2 <- function (d, sites, surveys = NULL, selection = NULL, traj.colors = NULL, traj.lwd = NULL,
                             axes = c(1, 2), survey.labels = FALSE, ...) {
  siteIDs = unique(sites)
  nsite = length(siteIDs)
  if (is.null(selection))
    selection = 1:nsite
  else {
    if (is.character(selection))
      selection = (siteIDs %in% selection)
  }
  selIDs = siteIDs[selection]
  xp = d[sites %in% selIDs, axes[1]]
  yp <- d[sites %in% selIDs, axes[2]]
  plot(xp, yp, type = "n", xlab = paste0("Axis ", axes[1]),
       ylab = paste0("Axis ", axes[2]), ...)
  sitesred = sites[sites %in% selIDs]
  if (!is.null(surveys))
    surveysred = surveys[sites %in% selIDs]
  else surveysred = NULL
  for (i in 1:length(selIDs)) {
    ind_surv = which(sitesred == selIDs[i])
    if (!is.null(surveysred))
      ind_surv = ind_surv[order(surveysred[sitesred ==
                                             selIDs[i]])]
    for (t in 1:(length(ind_surv) - 1)) {
      niini = ind_surv[t]
      nifin = ind_surv[t + 1]
      if (!is.null(traj.colors) & is.null(traj.lwd)){
        shape::Arrows(xp[niini], yp[niini], xp[nifin], yp[nifin],
                      col = traj.colors[i], arr.adj = 1, ...)
      } else if (!is.null(traj.colors) & !is.null(traj.lwd)){
        shape::Arrows(xp[niini], yp[niini], xp[nifin], yp[nifin],
                      col = traj.colors[i], arr.adj = 1, lwd = traj.lwd[i], ...)
      } else {shape::Arrows(xp[niini], yp[niini], xp[nifin], yp[nifin], arr.adj = 1)}

    }
  }
}

