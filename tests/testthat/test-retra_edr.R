test_that("returns an object of class 'RETRA'", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                         states = states, minSegs = 5)

  expect_s3_class(retra, "RETRA")
  expect_equal(attributes(retra$T1)$names,
               c("minSegs", "Segments", "Size", "Length",
                 "Link_distance", "Seg_density"))

  expect_type(retra$T1$minSegs, "double")
  expect_type(retra$T1$Segments, "character")
  expect_type(retra$T1$Size, "integer")
  expect_type(retra$T1$Length, "double")
  expect_s3_class(retra$T1$Link_distance, "data.frame")
  expect_type(retra$T1$Link_distance$Link, "character")
  expect_type(retra$T1$Link_distance$Distance, "double")
  expect_s3_class(retra$T1$Seg_density, "data.frame")
  expect_type(retra$T1$Seg_density$Density, "double")
  expect_type(retra$T1$Seg_density$kdTree_depth, "double")

})

test_that("returns same results for dSegs and coordSegs", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  # Segment dissimilarity
  dSegs <- ecotraj::segmentDistances(d = d, sites = trajectories, surveys = states)
  dSegs <- dSegs$Dseg
  seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", labels(dSegs))), "-")
  traj_Segs <- vapply(seg_components, function(iseg){
    iseg[1]
  }, character(1))
  state1_Segs <- vapply(seg_components, function(iseg){
    iseg[2]
  }, character(1))
  state2_Segs <- vapply(seg_components, function(iseg){
    iseg[3]
  }, character(1))

  # Segment coordinates
  ndim <- sum(table(trajectories) - 1) - 1
  set.seed(123)
  coordSegs <- smacof::mds(delta = dSegs, ndim = ndim, itmax = 300)
  coordSegs <- coordSegs$conf
  all(row.names(coordSegs) == labels(dSegs))

  # Compute RETRA-EDR
  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)
  retra_dSegs <- retra_edr(d = d, trajectories = trajectories,
                           states = states, minSegs = 5,
                           dSegs = dSegs, traj_Segs = traj_Segs,
                           state1_Segs = state1_Segs, state2_Segs = state2_Segs)
  retra_coordSegs <- retra_edr(d = d, trajectories = trajectories,
                               states = states, minSegs = 5,
                               dSegs = dSegs, coordSegs = coordSegs, traj_Segs = traj_Segs,
                               state1_Segs = state1_Segs, state2_Segs = state2_Segs)

  expect_equal(retra, retra_dSegs)
  expect_equal(retra, retra_coordSegs)
  expect_warning(retra_edr(d = d, trajectories = trajectories,
                           states = states, minSegs = 5,
                           dSegs = NULL, coordSegs = coordSegs, traj_Segs = traj_Segs,
                           state1_Segs = state1_Segs, state2_Segs = state2_Segs))

})

test_that("minSegs and eps works", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  nSegs <- sum(table(trajectories) - 1)
  retra_minSegs5_eps0 <- retra_edr(d = d, trajectories = trajectories, states = states,
                                   minSegs = 5)
  retra_minSegsMAX <- retra_edr(d = d, trajectories = trajectories, states = states,
                               minSegs = nSegs-1)
  retra_eps5 <- retra_edr(d = d, trajectories = trajectories, states = states,
                           minSegs = 5, eps = 5)

  expect_lte(length(retra_minSegsMAX), length(retra_minSegs5_eps0))
  expect_lte(length(retra_eps5), length(retra_minSegs5_eps0))

})

test_that("returns errors", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  nSegs <- sum(table(trajectories) - 1)

  # Segment dissimilarity
  dSegs <- ecotraj::segmentDistances(d = d, sites = trajectories, surveys = states)
  dSegs <- as.matrix(dSegs$Dseg)
  seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", labels(dSegs)[[1]])), "-")
  traj_Segs <- vapply(seg_components, function(iseg){
    iseg[1]
  }, character(1))
  state1_Segs <- vapply(seg_components, function(iseg){
    iseg[2]
  }, character(1))
  state2_Segs <- vapply(seg_components, function(iseg){
    iseg[3]
  }, character(1))

  # Segment coordinates
  ndim <- sum(table(trajectories) - 1) -1
  set.seed(123)
  coordSegs <- smacof::mds(delta = dSegs, ndim = ndim, itmax = 300)
  coordSegs <- coordSegs$conf
  all(row.names(coordSegs) == labels(dSegs)[[1]])

  expect_error(retra_edr(d = list(d), trajectories = trajectories, states = states, minSegs = 5),
               regexp = "'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  expect_error(retra_edr(d = d, trajectories = unique(trajectories), states = states, minSegs = 5),
               regexp = "The length of 'trajectories' must be equal to both dimensions in 'd'.")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = unique(states), minSegs = 5),
               regexp = "The length of 'states' must be equal to both dimensions in 'd'.")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, traj_Segs = traj_Segs),
               regexp = "To use 'dSegs' or 'coordSegs', you must provide values for 'traj_Segs', 'state1_Segs', and 'state2_Segs'.")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, traj_Segs = c("A", traj_Segs[2:length(traj_Segs)]),
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "Each value in 'traj_Segs' must be included in 'trajectories'.")
  expect_error(retra_edr(d = d, trajectories = c("A", trajectories[2:length(trajectories)]), states = states, minSegs = 5,
                         dSegs = dSegs, traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "Each value in 'trajectories' must be included in 'traj_Segs'.")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = c(100, state2_Segs[2:length(state2_Segs)])),
               regexp = "Each value in 'state1_Segs' and 'state2_Segs' must be included in 'states' for each corresponding site.")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = as.integer(c(100, states[2:length(states)])), minSegs = 5,
                         dSegs = dSegs, traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "Each value in 'states' must be included in at least one: 'state1_Segs' or 'state2_Segs', for each corresponding site.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = data.frame(dSegs), traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "'dSegs' must be a symmetric dissimilarity matrix or an object of class 'dist'.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs[2:(nrow(dSegs)), 2:(ncol(dSegs))], traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "The dimensions of 'dSegs' do not coincide with the total number of segments expected from 'trajectories'.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, traj_Segs = traj_Segs[2:length(traj_Segs)],
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "The length of 'traj_Segs', 'state1_Segs', and 'state2_Segs' must be equal to both dimensions in 'dSegs'.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, coordSegs = data.frame(coordSegs), traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "matrix containing the coordinates of trajectory segments")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, coordSegs = coordSegs[2:(nrow(coordSegs)), ], traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "The number of rows in 'coordSegs' do not coincide with the total number of segments expected from 'trajectories'.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, coordSegs = as.matrix(cbind(coordSegs, data.frame(X = 1:nrow(coordSegs)), Y = 1:nrow(coordSegs))), traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs),
               regexp = "The number of columns in 'coordSegs'")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                         dSegs = dSegs, coordSegs = coordSegs[, 1:5], traj_Segs = traj_Segs,
                         state1_Segs = state1_Segs, state2_Segs = state2_Segs, Dim = 6),
               regexp = "'Dim' cannot be greater than the number of columns in 'coordSegs'.")

  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5, Dim = 125),
               regexp = "'Dim' must be an integer in the range 1:nSegs")
  expect_error(retra_edr(d = d, trajectories = trajectories, states = states, minSegs = nSegs))


})

