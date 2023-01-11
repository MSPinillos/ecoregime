#' Representative trajectories in Ecological Dynamic Regimes (RETRA-EDR)
#'
#' @description
#' `retra_edr()` applies the algorithm RETRA-EDR (Sánchez-Pinillos et al.) to identify
#' representative trajectories summarizing the main dynamical patters of an Ecological
#' Dynamic Regime (EDR).
#'
#' @param d Either a symmetric matrix or an object of class [`dist`] containing the
#' dissimilarities between each pair of states of all trajectories in the EDR.
#' @param trajectories Vector indicating the trajectory or site to which each
#' state in `d` belongs.
#' @param states Vector of integers indicating the order of the states in `d` for
#' each trajectory.
#' @param minSegs Integer indicating the minimum number of segments in a region
#' of the EDR represented by a segment of the representative trajectory.
#' @param dSegs Either a symmetric matrix or an object of class [`dist`] containing
#' the dissimilarities between every pair of trajectory segments (see Details).
#' @param coordSegs Matrix containing the coordinates of trajectory segments (rows)
#' in each axis (columns) of an ordination space (see Details).
#' @param traj_Segs Vector indicating the trajectory to which each segment in `dSeg`
#' and/or `coordSegs` belongs. Only required if `dSegs` or `coordSegs` are not `NULL`.
#' @param state1_Segs Vector indicating the initial state of each segment in `dSegs`
#' and/or `coordSegs` according to the values given in `states`. Only required if
#' `dSegs` or `coordSegs` are not `NULL`.
#' @param state2_Segs Vector indicating the final state of each segment in `dSegs`
#' and/or `coordSegs` according to the values given in `states`. Only required if
#' `dSegs` or `coordSegs` are not `NULL`.
#' @param Dim Optional integer indicating the number of axes considered to
#' partition the segment space and generate a kd-tree. By default (`Dim = NULL`),
#' all axes are considered.
#' @param eps Numeric value indicating the minimum length in the axes of the segment
#' space to be partitioned when the kd-tree is generated. If `eps = 0` (default),
#' partitions are made regardless of the size.
#'
#' @details
#' The algorithm RETRA-EDR is based on a partition-and-group approach by which it
#' identifies regions densely crossed by ecological trajectories in an EDR, selects
#' a representative segment in each dense region, and joins the representative
#' segments by a set of artificial `Links` to generate a network of representative
#' trajectories. For that, RETRA-EDR splits the trajectories of the EDR into
#' segments and uses an ordination space generated from a matrix compiling the
#' dissimilarities between trajectory segments. Dense regions are identified by
#' applying a kd-tree to the ordination space.
#'
#' By default, RETRA-EDR calculates segment dissimilarities following the approach
#' by De Caceres et al. (2019) and applies metric multidimensional scaling (MDS,
#' Borg and Groenen, 2005) to generate the ordination space. It is possible to use
#' other dissimilarity metrics and/or ordination methods and reduce the computational
#' time by indicating the dissimilarity matrix and the coordinates of the segments
#' in the ordination space through the arguments `dSegs` and `coordSegs`, respectively.
#'
#' * If `!is.null(dSegs)` and `is.null(coordSegs)`, RETRA-EDR is computed by
#' applying MDS to `dSegs`.
#' * If `!is.null(dSegs)` and `!is.null(coordSegs)`, RETRA-EDR is directly computed
#' from the coordinates provided in `coordSegs` and representative segments are
#' identified using `dSegs`. `coordSegs` should be calculated by the user from
#' `dSegs`.
#' * If `is.null(dSegs)` and `!is.null(coordSegs)` (not recommended), RETRA-EDR
#' is directly computed from the coordinates provided in `coordSegs`. As `dSegs`
#' is not provided, `retra_edr()` assumes that the ordination space is metric and
#' identifies representative segments using the Euclidean distance.
#'
#' @return
#' The function `retra_edr()` returns an object of class `RETRA`, which is a list
#' of length equal to the number of representative trajectories identified. For
#' each trajectory, the following information is returned:
#' \describe{
#' \item{`minSegs`}{Value of the `minSegs` parameter.}
#' \item{`Segments`}{Vector of strings including the sequence of segments forming the
#' representative trajectory. Each segment is identified by a string of the form
#' `traj[st1-st2]`, where `traj` is the identifier of the original trajectory to
#' which the segment belongs and `st1` and `st2` are identifiers of the initial
#' and final states defining the segment.}
#'  \item{`Size`}{Numeric value indicating the number of states forming the representative
#' trajectory.}
#' \item{`Length`}{Numeric value indicating the length of the representative trajectory,
#' calculated as the sum of the dissimilarities in `d` between every pair of
#' consecutive states.}
#' \item{`Link_distance`}{Data frame of two columns indicating artificial links between
#' representative segments (`Link`) and the dissimilarity between the connected
#' states (`Distance`). When two representative segments are linked by a common
#' state or by two consecutive states of the same trajectory, the link distance
#' is zero or equal to the length of a real segment, respectively. In both cases,
#' the link is not considered in the returned data frame.}
#' \item{`Seg_density`}{Data frame of two columns and one row for each representative
#' segment. `Density` contains the number of segments in the leaf of the kd-tree
#' represented by each segment. `kdTree_depth` contains the depth of the kd-tree
#' for each leaf represented by the corresponding segment. That is, the number of
#' partitions of the ordination space until finding a region with `minSegs` segments
#' or less.}
#' }
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
#' @seealso
#' [`summary()`] for summarizing the characteristics of the set of representative
#' trajectories.
#'
#' [`plot()`] for plotting representative trajectories in an ordination space
#' representing the state space of the EDR.
#'
#' [`define_retra()`] for defining representative trajectories from a subset of
#' segments or trajectory features.
#'
#' @export
#'
#' @examples
#' # Example 1 -----------------------------------------------------------------
#' # Identify representative trajectories from state dissimilarities
#'
#' # Calculate state dissimilarities (Bray-Curtis) from species abundances
#' abundance <- data.frame(EDR_data$EDR1$abundance)
#' d <- vegan::vegdist(abundance[, -c(1:3)], method = "bray")
#'
#' # Identify the trajectory (or site) and states in d
#' trajectories <- abundance$traj
#' states <- as.integer(abundance$state)
#'
#' # Compute RETRA-EDR
#' RT1 <- retra_edr(d = d, trajectories = trajectories, states = states,
#'                  minSegs = 5)
#'
#' # Example 2 -----------------------------------------------------------------
#' # Identify representative trajectories from segment dissimilarities
#'
#' # Calculate segment dissimilarities using the Hausdorff distance
#' dSegs <- ecotraj:: segmentDistances(d = d, sites = trajectories,
#'                                     surveys = states,
#'                                     distance.type = "Hausdorff")
#' dSegs <- dSegs$Dseg
#'
#' # Identify the trajectory (or site) and states in dSegs:
#' # Split the labels of dSegs (traj[st1-st2]) into traj, st1, and st2
#' seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", labels(dSegs))), "-")
#' traj_Segs <- sapply(seg_components, "[", 1)
#' state1_Segs <- as.integer(sapply(seg_components, "[", 2))
#' state2_Segs <- as.integer(sapply(seg_components, "[", 3))
#'
#' # Compute RETRA-EDR
#' RT2 <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
#'                 dSegs = dSegs, traj_Segs = traj_Segs,
#'                 state1_Segs = state1_Segs, state2_Segs = state2_Segs)
#'
#'
retra_edr <- function (d, trajectories, states, minSegs,
                       dSegs = NULL, coordSegs = NULL,
                       traj_Segs = NULL, state1_Segs = NULL, state2_Segs = NULL,
                       Dim = NULL, eps = 0) {

  C1 = C3 = NULL # due to NSE notes in R CMD check

  ## WARNING MESSAGES ----------------------------------------------------------

  # Check the format for d, trajectories, and states
  if (all(!is.matrix(d), !dendextend::is.dist(d)) |
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

  # Check that traj_Segs and state_Segs are provided
  if (any(!is.null(dSegs), !is.null(coordSegs))) {

    if (any(is.null(traj_Segs), is.null(state1_Segs), is.null(state2_Segs))) {
      stop("To use 'dSegs' or 'coordSegs', you must provide values for 'traj_Segs', 'state1_Segs', and 'state2_Segs'.")
    }
    if (!all(traj_Segs %in% trajectories)) {
      stop("Each value in 'traj_Segs' must be included in 'trajectories'.")
    }
    if (!all(trajectories %in% traj_Segs)) {
      stop("Each value in 'trajectories' must be included in 'traj_Segs'.")
    }
    if (!all(unique(c(paste0(traj_Segs, state1_Segs), paste0(traj_Segs, state2_Segs))) %in% paste0(trajectories, states)) &
        length(unique(c(paste0(traj_Segs, state1_Segs), paste0(traj_Segs, state2_Segs)))) != length(paste0(trajectories, states))) {
      stop("Each value in 'state1_Segs' and 'state2_Segs' must be included in 'states' for each corresponding site.")
    }
    if (!all(paste0(trajectories, states) %in% unique(c(paste0(traj_Segs, state1_Segs), paste0(traj_Segs, state2_Segs))))) {
      stop("Each value in 'states' must be included in at least one: 'state1_Segs' or 'state2_Segs', for each corresponding site.")
    }
  }

  # Check the format of dSegs
  if (!is.null(dSegs)) {
    if (all(!is.matrix(dSegs), !dendextend::is.dist(dSegs)) |
        dim(dSegs)[1] != dim(dSegs)[2]) {
      stop("'dSegs' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
    }
    if (dim(dSegs)[1] != sum(table(trajectories) - 1)) {
      stop("The dimensions of 'dSegs' do not coincide with the total number of segments expected from 'trajectories'.")
    }
    if (any(nrow(as.matrix(dSegs)) != length(traj_Segs),
            nrow(as.matrix(dSegs)) != length(state1_Segs),
            nrow(as.matrix(dSegs)) != length(state2_Segs))) {
      stop("The length of 'traj_Segs', 'state1_Segs', and 'state2_Segs' must be equal to both dimensions in 'dSegs'.")
    }

    # Set names in dSegs
    dSegs <- as.matrix(dSegs)
    dimnames(dSegs) <- list(paste0(traj_Segs, "[", state1_Segs, "-", state2_Segs, ']'),
                            paste0(traj_Segs, "[", state1_Segs, "-", state2_Segs, ']'))

  }

  # Check the format of coordSegs
  if (!is.null(coordSegs)) {
    if(!is.matrix(coordSegs)) {
      stop("'coordSegs' must be a matrix containing the coordinates of trajectory segments (rows) in each axis (columns) of an ordination space.")
    }
    if (nrow(coordSegs) != sum(table(trajectories) - 1)) {
      stop("The number of rows in 'coordSegs' do not coincide with the total number of segments expected from 'trajectories'.")
    }
    if (ncol(coordSegs) > nrow(coordSegs)) {
      stop("The number of columns in 'coordSegs' (axes of the ordination space) cannot be greater than the number of rows (trajectory segments)")
    }
    if (any(nrow(coordSegs) != length(traj_Segs),
            nrow(coordSegs) != length(state1_Segs),
            nrow(coordSegs) != length(state2_Segs))) {
      stop("The length of 'traj_Segs', 'state1_Segs', and 'state2_Segs' must be equal to the number of rows in 'coordSegs'.")
    }
    rownames(coordSegs) <- paste0(traj_Segs, "[", state1_Segs, "-", state2_Segs, ']')

    if (!is.null(Dim)) {
      if (Dim > ncol(coordSegs)) {
        stop("'Dim' cannot be greater than the number of columns in 'coordSegs'.")
      }
    }
    if (is.null(Dim)){
      Dim <- ncol(coordSegs)
    }
  }

  # Check the format and value of minSegs and Dim
  if (minSegs %% 1 != 0 | !(minSegs %in% 1:(sum(table(trajectories) - 1)-1))) {
    stop("'minSegs' must be an integer in the range 1:(nSegs-1) (nSegs = total number of trajectory segments).")
  }
  if (!is.null(Dim)) {
    if (Dim %% 1 != 0 | !(Dim %in% 1:sum(table(trajectories) - 1))) {
      stop("'Dim' must be an integer in the range 1:nSegs (nSegs = total number of trajectory segments).")
    }
  }

  ## SEGMENT DISSIMILARITY AND SEGMENT SPACE -----------------------------------

  # Calculate segment distances
  if (all(is.null(dSegs), is.null(coordSegs))){
    dSegs <- ecotraj::segmentDistances(d = d, sites = trajectories, surveys = states)
    dSegs <- dSegs$Dseg
  }

  # Calculate the coordinates in a MDS
  if (is.null(Dim)) {
    Dim <- sum(table(trajectories) - 1)
  }
  if (Dim == sum(table(trajectories) - 1)) {
    ndim <- Dim - 1
  } else {
    ndim <- Dim
  }

  if (is.null(coordSegs)) {
    set.seed(123)
    coordSegs <- smacof::mds(delta = dSegs, ndim = ndim, itmax = 300)
    coordSegs <- coordSegs$conf
  }

  # Store data in a list and set keys
  seg_tree <- list(`0` = data.table::data.table(coordSegs, keep.rownames = T))
  data.table::setkey(seg_tree[[1]])
  nSegs = nrow(seg_tree[[1]])
  if(Dim > 1e5) {
    warning("Dim > 1e5, the first 1e5 dimensions will be used")
  }
  Dim_col = rep(1:Dim, 1e5) + 1

  ## COMPUTE KD-TREE BASED ON MIDPOINTS ----------------------------------------

  # First partition
  if (nSegs > minSegs) {
    i = 1                                 # i = tree depth
    iSplit <- 1                           # First partition
    icol <- Dim_col[i]                    # col id corresponding to the partition axis
    iRange <- range(seg_tree[[iSplit]][[icol]])   # Min-max values of the partition axis

    # Calculate the threshold to split the data
    cutpoint <- mean(iRange)

    # Define the set to split
    toSplit <- seg_tree[[iSplit]]

    # Split toSplit and store in a list
    curSplit <- list(toSplit[toSplit[[icol]] <= cutpoint],
                     toSplit[toSplit[[icol]] > cutpoint])
    names(curSplit) <- paste0(names(seg_tree), 0:1)

    # Append curSplit to the root
    seg_tree <- append(seg_tree, curSplit)

    # Update nSegs, iRange, and toSplit for the next partition
    nSegs <- sapply(curSplit, nrow)
    iRange <- sapply(curSplit, function(jcurSplit) {
      range(jcurSplit[[Dim_col[(i+1)]]])
    })
    if (any(abs(diff(iRange)) > eps)) {
      toSplit <- curSplit[which(nSegs > minSegs & abs(diff(iRange)) > eps)]
    } else {
      toSplit <- integer()
    }
    len_toSplit <- length(toSplit)
  }

  # Next partitions
  while (len_toSplit > 0) {
    i = i+1
    icol = Dim_col[i]  # the axis used at sequential partitions changes in 1:Dim repeatedly

    for (iSplit in 1:len_toSplit) {

      # Define cut point
      cutpoint <- mean(iRange[, iSplit])

      # Split toSplit and store in a list
      curSplit <- list(toSplit[[iSplit]][toSplit[[iSplit]][[icol]] <= cutpoint],
                       toSplit[[iSplit]][toSplit[[iSplit]][[icol]] > cutpoint])
      names(curSplit) <- paste0(names(toSplit)[iSplit], 0:1)

      # Append curSplit the tree
      seg_tree <- append(seg_tree, curSplit)

      # Update nSegs and toSplit for the next partition
      nSegs <- sapply(curSplit, nrow)
      curRange <- sapply(curSplit, function(jcurSplit){range(jcurSplit[[Dim_col[(i+1)]]])})

      cur_toSplit <- which(nSegs > minSegs & abs(diff(curRange)) > eps)
      toSplit <- append(toSplit, curSplit[cur_toSplit])
      iRange <- cbind(iRange, curRange[, cur_toSplit])
    }

    # Remove data for nodes that are already split
    toSplit <- toSplit[-c(1:len_toSplit)]
    iRange <- as.matrix(iRange[, -c(1:len_toSplit)])

    # Update len_toSplit
    len_toSplit <- length(toSplit)

  }

  ## SELECT DENSE PARTITIONS (LEAVES) ------------------------------------------

  # Select partitions with > minSegs
  seg_tree <- seg_tree[lapply(seg_tree, nrow) > minSegs]
  segtree_nms <- names(seg_tree)
  leaves <- list()

  # Remove redundant data and keep the leaves of the tree
  while(length(segtree_nms) > 0){
    len_segtree <- length(segtree_nms)
    leaves <- append(leaves, seg_tree[segtree_nms[len_segtree]])
    redundant <- sapply(2:nchar(segtree_nms[len_segtree]), function(i){
      stringr::str_sub(segtree_nms[len_segtree], 1, -i)
    })

    segtree_nms <- segtree_nms[-which(segtree_nms %in% c(redundant, segtree_nms[len_segtree]))]
  }

  ## IDENTIFY MEDOIDS ----------------------------------------------------------

  leaf_dSegs <- list()
  id_medoid <- list()
  medoids <- list()

  if (is.null(dSegs)) {
    warning(cat("Representative trajectories were identified using Euclidean distances in 'coordSegs'. \nTo use a specific metric, provide a value for 'dSegs'."))
  }
  for (ipartition in seq_along(leaves)) {
    # Distance between all segments in each leaf
    if (is.null(dSegs)) {
      leaf_dSegs[[ipartition]] <- dist(leaves[[ipartition]][, 2:ncol(leaves[[ipartition]])])
    } else {
      leaf_dSegs[[ipartition]] <- as.dist(as.matrix(dSegs)[leaves[[ipartition]]$rn, leaves[[ipartition]]$rn])
    }


    # ID of the medoid segment in each leaf
    id_medoid[[ipartition]] <- GDAtools::medoids(D = leaf_dSegs[[ipartition]],
                                                 cl = rep(1, nrow(leaves[[ipartition]])))

    # Extract the information of medoids
    medoids[[ipartition]] <- leaves[[ipartition]][id_medoid[[ipartition]], ]

    # Identify tree node
    medoids[[ipartition]]$Node <- names(leaves)[ipartition]

  }

  medoids <- data.table::rbindlist(medoids)
  medoids_nms <- medoids[[1]]

  ## JOIN MEDOID SEGMENTS #---------------------------------------------------

  if (length(medoids_nms) <= 1) {
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", medoids_nms)), "-")
    itraj_states <- unlist(lapply(seg_components, function(iseg){
      c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
    }))
    repr_traj <- list(T1 = list(minSegs = minSegs,
                                Segments = medoids_nms,
                                Size = length(unique(itraj_states)),
                                Length = traj_length(traj_segs = medoids_nms, dState = d,
                                                     trajectories = trajectories, states = states),
                                Link_distance = NA,
                                Seg_density = data.frame(
                                  Density = vapply(medoids_nms, function (iseg) {
                                    node = medoids[which(medoids$rn == iseg), ]$Node
                                    density = nrow(leaves[[node]])
                                  }, numeric(1)),
                                  kdTree_depth = vapply(medoids_nms, function (iseg) {
                                    node = medoids[which(medoids$rn == iseg), ]$Node
                                    depth = nchar(node)-1
                                  }, numeric(1)))
    ))
  } else if (length(medoids_nms) > 1) {

    # State distances for medoids
    d <- as.matrix(d)
    medoid_sites <- sapply(medoids_nms, function(x) {
      unlist(strsplit(x, '\\['))[1]
    })
    idst <- which(trajectories %in% medoid_sites)
    d_medoids <- d[idst, idst]

    # Calculate segment dissimilarities for the medoids
    dSegs_medoids <- ecotraj::segmentDistances(d = d_medoids, sites = trajectories[idst], surveys = states[idst])
    dSegs_medoids <- as.matrix(dSegs_medoids$Dinifin)[medoids_nms, medoids_nms]

    # Identify the minimum distance between initial and final segment states
    nmeds <- length(medoids_nms)
    id_dist <- c("Dfinini", "Dinifin")
    linkD <- list(id_dist = matrix(0, nmeds, nmeds, dimnames = labels(dSegs_medoids)),
                  dSegs = matrix(0, nmeds, nmeds, dimnames = labels(dSegs_medoids)))
    for (r in 1:nmeds) {
      for (c in 1:nmeds) {
        linkD$id_dist[r, c] <- id_dist[which.min(c(dSegs_medoids[c, r], dSegs_medoids[r, c]))]
        linkD$dSegs[r, c] <- min(dSegs_medoids[r, c], dSegs_medoids[c, r])
      }
    }

    # Minimum spanning tree
    MST_medoids <- ape::mst(linkD$dSegs)

    # Pairs of linked segments
    sel_link <- which(MST_medoids == 1, arr.ind = T)
    links <- data.table::data.table(C1 = labels(MST_medoids)[[1]][sel_link[, 1]],
                                    C2 = labels(MST_medoids)[[2]][sel_link[, 2]],
                                    id_dist = linkD$id_dist[sel_link],
                                    dSegs = linkD$dSegs[sel_link])
    data.table::setkey(links)
    links <- links[id_dist == "Dfinini", c("C1", "C2")]

    # Construct representative trajectories
    Cp1 = 3
    temp <- merge(links, links, by.x = paste0("C", Cp1-1), by.y = paste0("C", Cp1-2))
    names(temp)[ncol(temp)] <- paste0("C", Cp1)
    temp <- temp[C1 != C3]
    data.table::setkey(temp)
    links <- rm_dup(links)

    while(nrow(temp) > 0 & ncol(temp) <= nmeds){
      # temp <- rm_dup(temp)
      temp <- temp[!duplicated(temp)]
      links <- merge(links, temp, all.x = T, by = data.table::.EACHI)
      data.table::setkey(links)
      Cp1 <- Cp1 + 1
      cols <- c(paste0("C", c(Cp1-2, Cp1-1)))
      temp <- merge(temp, temp[, cols, with = F],
                    by.x = paste0("C", Cp1-1), by.y = paste0("C", Cp1-2))
      names(temp)[ncol(temp)] <- paste0("C", Cp1)
      cols <- c(paste0("C", Cp1-2), paste0("C", Cp1))
      temp <- temp[temp[[cols[1]]] != temp[[cols[2]]]]
      data.table::setkey(temp)
      links <- rm_dup(links)
    }

    # Remove potential trajectories passing twice by the same segment
    rm.rows <- which(sapply(1:ncol(data.table::transpose(links)), function(i){
      any(duplicated(na.omit(data.table::transpose(links)[[i]])))
    }))
    if(length(rm.rows > 0)){
      links <- links[-rm.rows, ]
    }

    # Representative trajectories
    repr_traj <- data.table::transpose(links)

    repr_traj <- lapply(setNames(repr_traj, paste0("T", 1:ncol(repr_traj))), function(icol){
      ReprSegs = icol[!is.na(icol)]
    })

    repr_traj <- lapply(repr_traj, function (itraj) {
      seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", itraj)), "-")
      itraj_states <- unlist(lapply(seg_components, function(iseg){
        c(paste0(iseg[1], "_", iseg[2]), paste0(iseg[1], "_", iseg[3]))
      }))
      Link_distance <-  data.frame(
        Link = vapply(1:(length(itraj) - 1), function (ilink) {
          paste0(itraj[ilink:(ilink+1)], collapse = " - ")
        }, character(1)),
        Distance = vapply(1:(length(itraj) - 1), function (ilink) {
          linkD$dSegs[itraj[ilink], itraj[ilink + 1]]
        }, numeric(1)),
        Real = vapply(1:(length(itraj) - 1), function (ilink) {
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

      traj_attr <- list(minSegs = minSegs,
                        Segments = itraj,
                        Size = length(unique(itraj_states)),
                        Length = traj_length(traj_segs = itraj, dState = d,
                                             trajectories = trajectories, states = states),

                        Link_distance = Link_distance,
                        Seg_density = data.frame(
                          Density = vapply(itraj, function (iseg) {
                            node = medoids[which(medoids$rn == iseg), ]$Node
                            density = nrow(leaves[[node]])
                          }, numeric(1)),
                          kdTree_depth = vapply(itraj, function (iseg) {
                            node = medoids[which(medoids$rn == iseg), ]$Node
                            depth = nchar(node)-1
                          }, numeric(1)))
      )
      return(traj_attr)
    })

  }

  class(repr_traj) <- "RETRA"

  return(repr_traj)
}

# Auxiliary function to remove duplicate trajectories
rm_dup <- function(linksORtemp){
  Traj <- sapply(1:nrow(linksORtemp), function(r){
    paste0(unlist(na.omit(data.table::transpose(linksORtemp[r, ]))), collapse = "-")
  })
  rm.replicates <- integer()
  for(t in seq_along(Traj)){
    if(sum(grepl(Traj[t], Traj[-t], fixed = T) > 0)) {
      rm.replicates <- c(rm.replicates, t)
    }
  }
  if(length(rm.replicates) > 0){
    linksORtemp <- linksORtemp[-rm.replicates, ]
  }

  return(linksORtemp)
}


# Auxiliary function to calculate the length of the representative trajectory
traj_length <- function(traj_segs, dState, trajectories, states){
  traj_segs2 = gsub("\\]", "", gsub("\\[", "-", traj_segs))

  traj_st <- character()
  for (iSeg in strsplit(traj_segs2, "-")) {
    traj_st <- c(traj_st,
                 paste0(iSeg[1], "_", iSeg[2]),
                 paste0(iSeg[1], "_", iSeg[3]))
  }

  ind_d <- numeric()
  for (iSt in traj_st) {
    ind_d <- c(ind_d, which(paste0(trajectories, "_", states) %in% iSt))
  }

  i = 1
  length_RT <- numeric()
  while(i <= length(ind_d)-1){
    length_RT <- c(length_RT, as.matrix(dState)[ind_d[i], ind_d[i+1]])
    i = i+1
  }
  length_RT <- sum(length_RT)

  return(length_RT)

}


