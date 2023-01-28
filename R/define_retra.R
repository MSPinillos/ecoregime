#' Define representative trajectories from trajectory features
#'
#' @description
#' Generate an object of class `RETRA` from a data frame containing trajectory
#' states to define representative trajectories in Ecological Dynamic Regimes (EDR).
#'
#' @param data A data frame indicating identifiers for the new representative
#' trajectories, the individual trajectories or sites to which the states belong,
#' the order of the states in the individual trajectories, and the identifier of
#' the representative trajectory to which the states belong (only if `!is.null(retra)`).
#' Alternatively, 'data' can be a vector or a list of character vectors including
#' the sequence of segments. See Details for further
#' clarifications to define `data`.
#' @param d Either a symmetric matrix or an object of class [`dist`] containing the
#' dissimilarities between each pair of states of all trajectories in the EDR.
#' If `NULL` (default), the length (`Length`) of the representative trajectories and the
#' distances between states of different trajectories or sites (`Link_distance`)
#' are not calculated.
#' @param trajectories Only needed if `!is.null(d)`. Vector indicating the
#' trajectory or site to which each state in `d` belongs.
#' @param states Only needed if `!is.null(d)`. Vector of integers indicating the
#' order of the states in `d` for each trajectory.
#' @param retra Object of class `RETRA` returned from [retra_edr()]. If `NULL`
#' (default), `minSegs` and `Seg_density` are not provided.
#'
#' @details
#' Each representative trajectory returned by the function [retra_edr()] corresponds
#' to the longest sequence of representative segments that can be linked according
#' to the criteria defined in the RETRA-EDR algorithm (Sánchez-Pinillos et al.).
#' One could be interested in splitting the obtained trajectories, considering
#' only a fraction of the returned trajectories, or defining representative
#' trajectories following different criteria than those in RETRA-EDR.
#' The function `define_retra()` allows generating an object of class `RETRA` that
#' can be used in other functions of `ecoregime` (e.g., [plot()]).
#'
#' For that, it is necessary to provide information about the set of segments or
#' trajectory states that form the new representative trajectory through the
#' argument `data`:
#'
#' - `data` can be defined as a **data frame** with as many rows as the number of
#' states in all representative trajectories and the following columns:
#' \describe{
#' \item{`RT`}{String indicating the identifier of the new representative trajectories.
#' Each identifier needs to appear as many times as the number of states forming
#' each representative trajectory.}
#' \item{`RT_traj`}{Vector indicating the individual trajectory in the EDR to which
#' each state of the representative trajectory belongs.}
#' \item{`RT_states`}{Vector of integers indicating the identifier of the states
#' forming the representative trajectories. Each integer must refer to the order
#' of the states in the individual trajectories of the EDR to which they belong.}
#' \item{`RT_retra`}{Only if the new trajectories are defined from representative
#' trajectories returned by [retra_edr()] (when `!is.null(retra)`). Vector of strings
#' indicating the representative trajectory in `retra` to which each state belongs.}
#' }
#'
#' - Alternatively, `data` can be defined as a **vector** or a **list of character
#' vectors** containing the sequence of segments in one representative trajectory
#' (then, `data` is a vector of characters) or in a set of representative trajectories
#' (list of character vectors with as many elements as the number of representative
#' trajectories). In any case, each segment needs to be specified in the form
#' `traj[st1-st2]`, where `traj` is the identifier of the original trajectory to
#' which the segment belongs and `st1` and `st2` are identifiers of the initial
#' and final states defining the segment. If only one state of an individual
#' trajectory is considered to form the representative trajectory, the corresponding
#' segment needs to be defined as `traj[st-st]`.
#'
#' @return
#' An object of class `RETRA`, which is a list of length equal to the number of
#' representative trajectories defined. For each trajectory, the following
#' information is returned:
#' \describe{
#' \item{`minSegs`}{Value of the `minSegs` parameter used in [retra_edr()]. If `retra`
#' is `NULL`, `minSegs` = `NA`.}
#' \item{`Segments`}{Vector of strings including the sequence of segments forming the
#' representative trajectory. Each segment is identified by a string of the form
#' `traj[st1-st2]`, where `traj` is the identifier of the original trajectory to
#' which the segment belongs and `st1` and `st2` are identifiers of the initial
#' and final states defining the segment. The same format `traj[st1-st2]` is
#' maintained when only one state of an individual trajectory is considered
#' (`st1` = `st2`). `traj`, `st1`, and `st2` are recycled from `data`.}
#' \item{`Size`}{Integer indicating the number of states forming the representative
#' trajectory.}
#' \item{`Length`}{Numeric value indicating the length of the representative trajectory,
#' calculated as the sum of the dissimilarities in `d` between every pair of
#' consecutive states. If `d` is `NULL`, `Length` = `NA`.}
#' \item{`Link_distance`}{Data frame of two columns indicating artificial links between
#' two segments (`Link`) and the dissimilarity between the connected states
#' (`Distance`). When two representative segments are linked by a common state or
#' by two consecutive states of the same trajectory, the link distance is zero or
#' equal to the length of a real segment, respectively. In both cases, the link
#' is not considered in the returned data frame. If `d` is `NULL`, `Link_distance`
#' = `NA`.}
#' \item{`Seg_density`}{Data frame of two columns and one row for each representative
#' segment. `Density` contains the number of segments in the leaf of the kd-tree
#' represented by each segment. `kdTree_depth` contains the depth of the kd-tree
#' for each leaf represented by the corresponding segment. That is, the number of
#' partitions of the ordination space until finding a region with `minSegs` segments
#' or less. If `retra` is `NULL`, `Seg_density` = `NA`.}
#' }
#'
#' @author Martina Sánchez-Pinillos, CNRS, Univ. Montpellier
#'
#' @seealso
#' [retra_edr()] for identifying representative trajectories in EDRs through
#' RETRA-EDR.
#'
#' [`summary()`] for summarizing the characteristics of the set of representative
#' trajectories.
#'
#' [plot()] for plotting representative trajectories in an ordination space
#' representing the state space of the EDR.
#'
#' @export
#'
#' @examples
#' # Example 1 -----------------------------------------------------------------
#' # Define representative trajectories from the outputs of retra_edr().
#'
#' # Identify representative trajectories using retra_edr()
#' d <- EDR_data$EDR1$state_dissim
#' trajectories <- EDR_data$EDR1$abundance$traj
#' states <- EDR_data$EDR1$abundance$state
#' old_retra <- retra_edr(d = d, trajectories = trajectories, states = states,
#'                        minSegs = 5)
#'
#' # retra_edr() returns three representative trajectories
#' old_retra
#'
#' # Keep the last five segments of trajectories "T2" and "T3"
#' selected_segs <- old_retra$T2$Segments[4:length(old_retra$T2$Segments)]
#'
#' # Identify the individual trajectories for each state...
#' selected_segs   # each value represents traj[st1-st2]
#' selected_traj <- rep(c(15, 4, 4, 1, 14), each = 2) # This is "traj" in traj[st1-st2]
#'
#' # ...and the states (in the same order than the representative trajectory).
#' selected_states <- c(1, 2, 2, 3, 3, 4, 1, 2, 2, 3) # This is "st1, st2" in traj[st1-st2]
#'
#' # Generate the data frame with the format indicated in the documentation
#' df <- data.frame(RT = rep("A", length(selected_states)), # name of the new trajectory as "A"
#'                  RT_traj = selected_traj,
#'                  RT_states = as.integer(selected_states),
#'                  RT_retra = rep("T2", length(selected_states))) # it could be "T3" too
#'
#' # Remove duplicates (trajectory 4, state 3)
#' df <- unique(df)
#'
#' # Generate a RETRA object using define_retra()
#' new_retra <- define_retra(data = df,
#'                           d = d,
#'                           trajectories = trajectories,
#'                           states = states,
#'                           retra = old_retra)
#'
#' # Example 2 -----------------------------------------------------------------
#' # Define representative trajectories from sequences of segments
#'
#' # Select all segments in T1, split T2 into two new trajectories, and include
#' # a trajectory composed of states belonging to trajectories "5", "6", and "7"
#' data <- list(old_retra$T1$Segments,
#'              old_retra$T2$Segments[1:3],
#'              old_retra$T2$Segments[4:8],
#'              c("5[1-2]", "5[2-3]", "7[4-4]", "6[4-5]"))
#'
#' #' # Generate a RETRA object using define_retra()
#' new_retra <- define_retra(data = data,
#'                           d = d,
#'                           trajectories = trajectories,
#'                           states = states,
#'                           retra = old_retra)
#'
#' # Example 3 -----------------------------------------------------------------
#' # Define two representative trajectories from individual trajectories in EDR1.
#'
#' # Define trajectory "A" from states in trajectories 3 and 4
#' data_A <- data.frame(RT = rep("A", 4), # name of the new representative trajectory
#'                      RT_traj = c(3, 3, 4, 4), # identifier of the original trajectories
#'                      RT_states = c(1:2, 4:5)) # states in the original trajectories
#'
#' # Define trajectory "B" from states in trajectories 5, 6, and 7
#' data_B <- data.frame(RT = rep("B", 5), # name of the new representative trajectory
#'                      RT_traj = c(5, 5, 7, 6, 6), # identifier of the original trajectories
#'                      RT_states = c(1, 2, 4, 4, 5)) # states in the original trajectories
#'
#' # Compile data for both trajectories in a data frame
#' df <- rbind(data_A, data_B)
#' df$RT_states <- as.integer(df$RT_states)
#'
#' # Generate a RETRA object using define_retra()
#' new_retra <- define_retra(data = df, d = EDR_data$EDR1$state_dissim,
#'                           trajectories = EDR_data$EDR1$abundance$traj,
#'                           states = EDR_data$EDR1$abundance$state)
#'
#'
define_retra <- function(data, d = NULL, trajectories = NULL, states = NULL, retra = NULL) {

  ## CONVERT SEGMENTS INTO DATA FRAME ------------------------------------------

  # If data is a sequence of segments, generate a data.frame
  if (is.character(data)) {
    data <- segment_to_df(data, retra = retra)
  }
  if (is.list(data) & !is.data.frame(data)) {
    data.ls <- lapply(data, function(idata){
      data.df <- segment_to_df(segments = idata, retra = retra)
      data.df <- unique(data.df)
    })
    names(data.ls) <- vapply(data.ls, function(idata){
      if ("RT_retra" %in% names(idata)) {
        unique(idata$RT_retra)
      } else{
        unique(idata$RT)
      }
    }, character(1))
    for (id in unique(names(data.ls))) {
      len_idata <- length(data.ls[which(names(data.ls) == id)])
      if(len_idata > 1) {
        for (i in 1:len_idata) {
          data.ls[which(names(data.ls) == id)][i][[id]]$RT <- paste0(id, ".", i)
        }
      }
    }

    data <- data.frame(data.table::rbindlist(data.ls, fill = T))
  }

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for data
  if (sum(names(data) %in% c("RT", "RT_traj", "RT_states")) < 3) {
    stop("'data' must contain at least three columns named \"RT\", \"RT_traj\", and \"RT_states\".")
  }
  if (!is.integer(data$RT_states)) {
    stop("The column \"RT_states\" in 'data' needs to be of class integer.")
  }
  if (!is.null(retra)) {
    if (!("RT_retra" %in% names(data))) {
      warning("'data' does not contain a column named \"RT_retra\".")
      retra <- NULL
    }
  }

  # Check the format for d, trajectories, and states
  if (!is.null(d)) {
    if (all(!is.matrix(d), class(d) != "dist") |
        nrow(as.matrix(d)) != ncol(as.matrix(d))) {
      stop("'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
    }

    if (length(trajectories) != nrow(as.matrix(d))) {
      stop("The length of 'trajectories' must be equal to both dimensions in 'd'.")
    }
    if (length(states) != nrow(as.matrix(d))) {
      stop("The length of 'states' must be equal to both dimensions in 'd'.")
    }

  }

  ## ATTRIBUTE DEFINITION ------------------------------------------------------

  ID_RT <- unique(data$RT)
  nRT <- length(ID_RT)
  data$id <- 1:nrow(data)

  if (!is.null(retra)) {
    ID_RT_retra <- names(retra)
  }

  # Define attributes for each trajectory
  RT <- lapply(setNames(1:nRT, ID_RT), function(iT){
    ind_traj <- which(data$RT == ID_RT[iT])
    nStates <- length(ind_traj)
    idata <- data[ind_traj, ]

    # Some states can be unique in one trajectory. The segment will be defined
    # as traj[St-St]
    if (idata$RT_traj[1] != idata$RT_traj[2] |
        (idata$RT_traj[1] == idata$RT_traj[2] & idata$RT_states[2] - idata$RT_states[1] > 1)) {
      data <- rbind(data, data[ind_traj[1], ])
    }
    if (idata$RT_traj[nStates] != idata$RT_traj[nStates-1] |
        (idata$RT_traj[nStates] == idata$RT_traj[nStates-1] & idata$RT_states[nStates] - idata$RT_states[nStates-1] > 1)) {
      data <- rbind(data, data[ind_traj[nStates], ])
    }
    if (nStates > 2) {
      for (i in 2:(nStates-1)) {
        if (!(idata$RT_traj[i] %in% idata$RT_traj[c(i+1, i-1)]) |
            ((idata$RT_traj[i] == idata$RT_traj[(i-1)] & diff(idata$RT_states[c(i-1, i)]) > 1) |
             (idata$RT_traj[i] == idata$RT_traj[(i+1)] & diff(idata$RT_states[c(i, i+1)]) > 1))) {
          data <- rbind(data, data[ind_traj[i], ])
        }
      }
    }

    data <- data[sort(data$id), ]
    ind_traj <- which(data$RT == ID_RT[iT])
    nStates <- length(ind_traj)

    # Test that there is one RT_retra value
    if (!is.null(retra) && length(unique(data$RT_retra[ind_traj])) == 1) {
      id_retra <- unique(data$RT_retra[ind_traj])
    } else if (length(unique(data$RT_retra[ind_traj])) > 1) {
      warning("The column \"RT_retra\" in 'data' cannot refer to more than one trajectory in 'retra' for each trajectory specified in 'data$RT'.")
      id_retra <- NA
    } else if (is.null(retra)) {
      id_retra <- NA
    }

    # minSegs
    if (!is.null(retra) && id_retra %in% ID_RT_retra){
      minSegs <- retra[[id_retra]]$minSegs
    } else {
      minSegs <- NA
    }

    # Segments
    i <- 2
    Segments <- character()
    while (i <= nStates) {
      trajSt1 <- data$RT_traj[ind_traj][i-1]
      trajSt2 <- data$RT_traj[ind_traj][i]
      St1 <- data$RT_states[ind_traj][i-1]
      St2 <- data$RT_states[ind_traj][i]

      if (trajSt1 == trajSt2 && St2 - St1 <= 1) {
        Segments <- c(Segments, paste0(trajSt1, "[", St1, "-", St2, "]"))
      }
      i <- i+1
    }

    # Size
    Size <- length(unique(paste0(data$RT_traj[ind_traj], "_", data$RT_states[ind_traj])))

    # Length
    if (!is.null(d)) {
      Length <- traj_length(traj_segs = Segments, dState = d,
                            trajectories = trajectories,
                            states = states)
    } else if (is.null(d)) {
      Length <-  NA
    }

    # Link distance
    if (!is.null(d) & length(Segments) > 1) {
      traj_states <-paste0(trajectories, "_", states)
      seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", Segments)), "-")
      Link_distance = data.frame(
        Link = vapply(1:(length(Segments) - 1), function (ilink) {
          paste0(Segments[ilink:(ilink+1)], collapse = " - ")
        }, character(1)),

        Distance = vapply(1:(length(Segments) - 1), function (ilink) {
          iSt1 <- paste0(seg_components[[ilink]][1], "_", seg_components[[ilink]][3])
          iSt2 <- paste0(seg_components[[ilink+1]][1], "_", seg_components[[ilink+1]][2])
          if (iSt1 %in% traj_states & iSt2 %in% traj_states) {
            distance <- as.matrix(d)[which(traj_states == iSt1),
                                     which(traj_states == iSt2)]
          } else {
            distance <- NA
          }
          return(distance)
        }, numeric(1)),

        Real = vapply(1:(length(Segments) - 1), function (ilink) {
          traj1 <- seg_components[[ilink]][1]
          traj2 <- seg_components[[ilink+1]][1]
          iSt1 <- as.integer(seg_components[[ilink]][3])
          iSt2 <- as.integer(seg_components[[ilink+1]][2])
          if (traj1 == traj2 & diff(c(iSt1, iSt2)) <= 1) {T} else {F}
        }, logical(1))
      )
      Link_distance <- Link_distance[which(Link_distance$Real == F), c("Link", "Distance")]
      if (nrow(Link_distance) == 0){
        Link_distance <- NA
      } else {row.names(Link_distance) <- 1:nrow(Link_distance)}
    } else {
      Link_distance <- NA
    }

    # Seg density
    if (!is.null(retra) && id_retra %in% ID_RT_retra){
      irows <- which(rownames(retra[[id_retra]]$Seg_density) %in% Segments)
      Seg_density <- retra[[id_retra]]$Seg_density[irows, ]
    } else {
      Seg_density <- NA
    }

    # Compile all elements in a list
    new_retra <- list(minSegs = minSegs,
                      Segments = Segments,
                      Size = Size,
                      Length = Length,
                      Link_distance = Link_distance,
                      Seg_density = Seg_density)

    return(new_retra)
  })

  class(RT) <- "RETRA"
  return(RT)
}

# Auxiliary function to convert a vector of segments into a data.frame
segment_to_df <- function(segments, retra = retra) {
  seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", segments)), "-")
  RT_traj <- c(vapply(seg_components, function(iseg){
    c(iseg[1], iseg[1])
  }, character(2)))
  RT_states <- c(vapply(seg_components, function(iseg){
    c(as.integer(iseg[2]), as.integer(iseg[3]))
  }, integer(2)))

  if (!is.null(retra)) {
    retra_segs <- lapply(retra, "[", "Segments")
    ind_retra <- which(vapply(retra_segs, function(iretra){
      all(segments %in% iretra$Segments)
    }, logical(1)))
    if (length(ind_retra) == 0) {
      RT <- rep("newT", length(RT_states))
      data <- data.frame(RT = RT, RT_traj = RT_traj, RT_states = RT_states)
    }
    if (length(ind_retra) >= 1) {
      ind_retra <- ind_retra[1]
      RT_retra <- rep(names(retra)[ind_retra], length(RT_states))
      RT <- paste0(RT_retra, ".1")
      data <- data.frame(RT = RT, RT_traj = RT_traj, RT_states = RT_states,
                         RT_retra = RT_retra)
    }
  } else {
    RT <- rep("newT", length(RT_states))
    data <- data.frame(RT = RT, RT_traj = RT_traj, RT_states = RT_states)
  }

  data$dup <- sapply(1:(nrow(data)), function(irow){
    ifelse(any(duplicated(data[irow:(irow+1), ])), 1, 0)
  })

  data <- data[which(data$dup == 0), ]

  return(data)
}


