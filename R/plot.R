#' Plot Ecological Dynamic Regimes and representative trajectories
#'
#' @description
#' Plot individual trajectories belonging to an EDR in the state space and highlight
#' a set of representative trajectories, distinguishing representative segments
#' belonging to existent trajectories of the EDR and the links between consecutive
#' segments.
#'
#' @param x Object of class `RETRA` returned from functions [`retra_edr()`] or
#' [`define_retra()`].
#' @param d Symmetric matrix or `dist` object containing the dissimilarities
#' between each pair of states of all trajectories or data frame containing the
#' coordinates of all trajectory states in an ordination space.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `d` belongs.
#' @param states Vector of integers indicating the order of the states in `d` for
#' each trajectory.
#' @param select_RT Optional string indicating the name of representative trajectories
#' selected to be highlighted in the plot. Defaults NULL.
#' @param traj.colors A specification for the color of individual trajectories
#' (defaults "grey") or a vector with length equal to the number of trajectories
#' indicating the color for each one.
#' @param RT.colors A specification for the color of representative trajectories
#' (defaults "black").
#' @param sel.color A specification for the color of the selected representative
#' trajectory (defaults "red").
#' @param link.color A specification for the color of the links between trajectory
#' segments forming representative trajectories. Defaults the same color than
#' `RT.colors`.
#' @param link.lty The line type of the links between trajectory segments forming
#' representative trajectories. Defaults 2 = dashed (See [graphics::par]).
#' @param axes An integer vector indicating the pair of axes in the ordination
#' space to be plotted.
#' @param ... Arguments for function [shape::Arrows].
#'
#' @return
#' The function `plot.RETRA()` plots a set of individual trajectories and the
#' representative trajectories in an ordination space defined through `d` or
#' calculated by applying metric multidimensional scaling (mMDS) to `d`.
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
#' @note
#' This function uses a modified version of the trajectoryPlot() function in 'ecotraj'
#' (v.0.0.3; De Cáceres et al. 2019). The modification was done the 2022-12-18 and
#' allows using the function Arrows() in the 'shape' package to plot ecological
#' trajectories instead the function arrows() in base.
#'
#' @seealso
#' [`retra_edr()`] for identifying representative trajectories in EDRs applying
#' RETRA-EDR.
#' [`define_retra()`] for generating an object of class `RETRA` from a sequence
#' of trajectory states.
#' [`summary()`] for summarizing representative trajectories in EDRs.
#'
#' @export
#'
#' @examples
#' # Example 1 -----------------------------------------------------------------
#'
#' # `d` contains the dissimilarities between trajectory states
#' d = EDR_data$EDR1$state_dissim
#'
#' # `trajectories` and ` states` are defined according to `d` entries.
#' trajectories = EDR_data$EDR1$abundance$traj
#' states = EDR_data$EDR1$abundance$state
#'
#' # `x` defined from `retra_edr()`. We obtain three representative trajectories.
#' # We will select "T2" to be plotted in a different color.
#' RT <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)
#' summary(RT)
#'
#' # Plot individual trajectories in blue and representative trajectories in orange,
#' # except "T2", which will be displayed in green Artificial links will be displayed
#' # with a dotted line.
#' plot(x = RT, d = d, trajectories = trajectories, states = states, select_RT = "T2",
#'      traj.colors = "lightblue", RT.colors = "orange", sel.color = "darkgreen",
#'      link.lty = 3, main = "Representative trajectories in EDR1")
#'
#' # Example 2 -----------------------------------------------------------------
#'
#' # `d` contains the coordinates in an ordination space. For example, we use
#' # the coordinates of the trajectory states after applying a principal component
#' # analysis (PCA) to an abundance matrix.
#' abun <- EDR_data$EDR1$abundance
#' pca <- prcomp(abun[, -c(1:3)])
#' coord = data.frame(pca$x)
#'
#' # `trajectories` and ` states` are defined according to the abundance matrix
#' # used in the PCA
#' trajectories = EDR_data$EDR1$abundance$traj
#' states = EDR_data$EDR1$abundance$state
#'
#' # Instead of using the representative trajectories obtained from `retra_edr()`,
#' # we will define the set of "representative trajectories". For example, we can
#' # select the trajectories whose initial and final states are in the extremes
#' # of the first axis.
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
#' # Plot our selection of representative trajectories with the default graphic
#' # values
#' plot(x = RT_retra, d = coord, trajectories = trajectories, states = states,
#'      main = "Extreme trajectories in EDR1")
#'

plot.RETRA <- function (x, d, trajectories, states, select_RT = NULL,
                        traj.colors = NULL, RT.colors = NULL, sel.color = NULL,
                        link.color = NULL, link.lty = 2, axes = c(1, 2), ...) {

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the formats
  if (!is.data.frame(d)) {
    if (all(!is.matrix(d), !dendextend::is.dist(d)) |
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

  ## RETRESENTATIVE TRAJECTORIES -----------------------------------------------

  if (is(x, "RETRA")) {
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
        RT_states <- c(RT_states, iseg[2:3])
      }
    }

    if(sum(!paste0(RT_traj, RT_states) %in% paste0(trajectories, states)) > 0){
      stop("The states in forming representative trajectories (identified through 'Segments' in x) need to be included in 'trajectories' and 'states'.")
    }
  }

  if (!is.null(select_RT) && !(select_RT %in% RT_names)) {
    stop(cat("select_RT = \"", select_RT, "\"; \"", select_RT, "\" needs to be included in 'x' or 'x$RT_names'", sep = ""))
  }


  # STATE COORDINATES ----------------------------------------------------------

    # Coordinates MDS
  if (dendextend::is.dist(d) || isSymmetric(as.matrix(d))) {
    warning(cat("Representative trajectories will be displayed in an ordination space generated through multidimensional scaling (MDS). You can avoid this step by providing state coordinates in the 'd' argument.", "\n"))
    statesMDS <- data.frame(smacof::mds(delta = d, ndim = ncol(d)-1,
                                        itmax = 300, verbose = T)$conf)
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
  if (nRT > 1) {
    # Place the selected RT to be displayed at the end
    ID_RT <- names(sort(table(RT_names)))
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

  # Identify links
  links <- lapply(1:nRT, function(iT){
    ind_traj <- which(RT_data$RT == ID_RT[iT])

    # False links
    nlinks <- nrow(RT_data[ind_traj, ]) - 1
    nlini <- RT_traj[ind_traj][seq(2, nlinks, by = 2)]
    nlfin <- RT_traj[ind_traj][seq(3, nlinks, by = 2)]
    ind_link <- which(nlini != nlfin)

    if (length(ind_link) > 0) {
      links <- data.frame(ID = RT_data$ID[ind_traj][c(2*ind_link, 2*ind_link+1)],
                          RT = ID_RT[iT],
                          Link = seq_along(ind_link),
                          Link_state = rep(1:2, each = length(ind_link)), row.names = NULL)
      links <- links[order(links$Link, links$Link_state), ]
    } else {links <- data.frame()}

    return(links)
  })
  links <- data.frame(data.table::rbindlist(links))

  # Coordinates of link states
  if (nrow(links) > 0) {
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
  plot(xp, yp, type = "n", asp = 1, xlab = paste0("Axis ", axes[1]), ylab = paste0("Axis ", axes[2]), ...)
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

# #' @export
# to_RETRAdf <- function(x){
#   if (!is.data.frame(x)) {
#     stop("'x' needs to be a data frame.")
#   }
#   if (sum(names(x) %in% c("RT_names", "RT_traj", "RT_states")) != 3) {
#     stop("'x' needs to include data for 'RT_names', 'RT_traj', 'RT_states'.")
#   }
#
#   class(x) <- "RETRAdf"
#   return(x)
# }