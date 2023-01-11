test_that("returns an object of class 'RETRA'", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  data <- data.frame(RT = c(rep("A", 6), rep("B", 2)),
                     RT_traj = c(1, 1, 1, 2, 3, 3, 4, 4),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5, 1, 2)))
  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  expect_s3_class(new_retra, "RETRA")
  expect_equal(attributes(new_retra$A)$names,
               c("minSegs", "Segments", "Size", "Length", "Link_distance", "Seg_density"))

  expect_type(new_retra$A$Segments, "character")
  expect_type(new_retra$A$Size, "integer")
  expect_type(new_retra$A$Length, "double")
  expect_s3_class(new_retra$A$Link_distance, "data.frame")
  expect_equal(new_retra$B$Link_distance, NA)
  expect_type(new_retra$A$Link_distance$Link, "character")
  expect_type(new_retra$A$Link_distance$Distance, "double")

})

test_that("returns an object of class 'RETRA' (data is a list)", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  data <- data.frame(RT = c(rep("newT.1", 6), rep("newT.2", 2)),
                     RT_traj = c(1, 1, 1, 2, 3, 3, 4, 4),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5, 1, 2)))
  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  data2 <- list(c("1[1-2]", "1[2-3]", "2[3-3]", "3[4-5]"),
               c("4[1-2]"))
  new_retra2 <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  expect_s3_class(new_retra2, "RETRA")
  expect_equal(attributes(new_retra2$newT.1)$names,
               c("minSegs", "Segments", "Size", "Length", "Link_distance", "Seg_density"))

  expect_type(new_retra2$newT.1$Segments, "character")
  expect_type(new_retra2$newT.1$Size, "integer")
  expect_type(new_retra2$newT.1$Length, "double")
  expect_s3_class(new_retra2$newT.1$Link_distance, "data.frame")
  expect_equal(new_retra2$newT.2$Link_distance, NA)
  expect_type(new_retra2$newT.1$Link_distance$Link, "character")
  expect_type(new_retra2$newT.1$Link_distance$Distance, "double")

  expect_equal(new_retra, new_retra2)

})

test_that("returns same results when data is defined from 'retra'", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  old_retra <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)

  # Generate data to generate identical RETRA
  data.ls <- lapply(seq_along(old_retra), function(itraj){
    Segments <- old_retra[[itraj]]$Segments
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", Segments)), "-")
    RT_traj <- vapply(seg_components, function(iseg){
      c(iseg[1], iseg[1])
    }, character(2))
    RT_traj <- c(RT_traj)
    RT_states <- vapply(seg_components, function(iseg){
      c(as.integer(iseg[2]), as.integer(iseg[3]))
    }, integer(2))
    RT_states <- c(RT_states)

    data <- data.frame(RT_traj = RT_traj, RT_states = RT_states)
    data$RT <- names(old_retra)[itraj]
    data$RT_retra <- names(old_retra)[itraj]

    return(unique(data))
  })

  data1 <- do.call(rbind, data.ls)
  new_retra1 <- define_retra(data = data1, d = d,
                             trajectories = trajectories, states = states,
                             retra = old_retra)

  data2 <- sapply(old_retra, '[', 'Segments')
  new_retra2 <- define_retra(data = data2, d = d,
                             trajectories = trajectories, states = states,
                             retra = old_retra)
  names(new_retra2) <- names(old_retra)

  expect_equal(new_retra1, new_retra2)

  # The deinition of the segments may change: RETRA-EDR only returns the representative
  # segments (e.g., T[1-2], T[3-4]); define_retra returns all possible segments
  # (T[1-2], T[2-3], T[3-4])
  old_retra_noSegs <- lapply(old_retra, "[", c("minSegs", "Size", "Length",
                                               "Link_distance", "Seg_density"))
  new_retra_noSegs <- lapply(new_retra1, "[", c("minSegs", "Size", "Length",
                                               "Link_distance", "Seg_density"))

  expect_equal(old_retra_noSegs, new_retra_noSegs)

})

test_that("attributes change for a selection of states", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  old_retra <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)

  data.ls <- lapply(seq_along(old_retra), function(itraj){
    Segments <- c(old_retra[[itraj]]$Segments[1],
                  old_retra[[itraj]]$Segments[length(old_retra[[itraj]]$Segments)])
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", Segments)), "-")
    RT_traj <- vapply(seg_components, function(iseg){
      c(iseg[1], iseg[1])
    }, character(2))
    RT_traj <- c(RT_traj)
    RT_states <- vapply(seg_components, function(iseg){
      c(as.integer(iseg[2]), as.integer(iseg[3]))
    }, integer(2))
    RT_states <- c(RT_states)

    data <- data.frame(RT_traj = RT_traj, RT_states)
    data$RT <- names(old_retra)[itraj]
    data$RT_retra <- names(old_retra)[itraj]

    return(unique(data))
  })

  data <- do.call(rbind, data.ls)
  new_retra <- define_retra(data = data, d = as.matrix(d),
                             trajectories = trajectories, states = states,
                             retra = old_retra)

  data2 <-  lapply(seq_along(old_retra), function(itraj){
    Segments <- c(old_retra[[itraj]]$Segments[1],
                  old_retra[[itraj]]$Segments[length(old_retra[[itraj]]$Segments)])
  })
  new_retra2 <- define_retra(data = data2, d = as.matrix(d),
                             trajectories = trajectories, states = states,
                             retra = old_retra)
  names(new_retra2) <- names(new_retra)
  expect_equal(new_retra, new_retra2)


  size_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Size - new_retra[[iRT]]$Size
  })
  length_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Length - new_retra[[iRT]]$Length
  })
  link_diff <- lapply(seq_along(old_retra), function(iRT){
    if (is.data.frame(new_retra[[iRT]]$Link_distance) & is.data.frame(old_retra[[iRT]]$Link_distance)) {
      max(new_retra[[iRT]]$Link_distance$Distance) - max(old_retra[[iRT]]$Link_distance$Distance)
    } else {0}
  })

  expect_true(all(size_diff >= 0))
  expect_true(all(length_diff >= 0))
  expect_true(all(link_diff >= 0))

})

test_that("defines a segment", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  dimnames(d) <- list(paste0(trajectories, "_", states),
                      paste0(trajectories, "_", states))

  data <- data.frame(RT = rep("newT", 2),
                     RT_traj = rep(28, 2),
                     RT_states = as.integer(c(1, 2)))
  retra <- define_retra(data = data, d = d,
                        trajectories = trajectories, states = states)
  data2 <- "28[1-2]"
  retra2 <- define_retra(data = data2, d = d,
                        trajectories = trajectories, states = states)

  expect_equal(retra, retra2)
  expect_equal(length(retra$newT$Segments), 1)
  expect_equal(retra$newT$Size, 2)
  expect_equal(retra$newT$Length, d["28_1", "28_2"])
  expect_equal(retra$newT$Link_distance, NA)

})

test_that("defines a complete trajectory", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  dimnames(d) <- list(paste0(trajectories, "_", states),
                      paste0(trajectories, "_", states))

  data <- data.frame(RT = rep("newT", 5),
                     RT_traj = rep(1, 5),
                     RT_states = as.integer(1:5))
  retra <- define_retra(data = data, d = d,
                        trajectories = trajectories, states = states)

  data2 <- c("1[1-2]", "1[2-3]", "1[3-4]", "1[4-5]")
  retra2 <- define_retra(data = data2, d = d,
                         trajectories = trajectories, states = states)

  expect_equal(retra, retra2)

  expect_equal(length(retra$newT$Segments), nrow(data)-1)
  expect_equal(retra$newT$Size, nrow(data))
  expect_equal(retra$newT$Length, sum(d["1_1", "1_2"], d["1_2", "1_3"],
                                      d["1_3", "1_4"], d["1_4", "1_5"]))
  expect_equal(retra$newT$Link_distance, NA)

})

test_that("defines a segment composed of states from two trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  dimnames(d) <- list(paste0(trajectories, "_", states),
                      paste0(trajectories, "_", states))

  data <- data.frame(RT = rep("newT", 2),
                     RT_traj = c(28, 30),
                     RT_states = as.integer(c(1, 2)))
  retra <- define_retra(data = data, d = d,
                        trajectories = trajectories, states = states)

  data2 <- c("28[1-1]", "30[2-2]")
  retra2 <- define_retra(data = data2, d = d,
                        trajectories = trajectories, states = states)

  expect_equal(retra, retra2)

  expect_equal(length(retra$newT$Segments), 2)
  expect_equal(retra$newT$Size, 2)
  expect_equal(retra$newT$Length, d["28_1", "30_2"])
  expect_equal(retra$newT$Link_distance$Distance, retra$newT$Length)

})

test_that("defines a sequence of states from different trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  dimnames(d) <- list(paste0(trajectories, "_", states),
                      paste0(trajectories, "_", states))

  data <- data.frame(RT = rep("newT", 4),
                     RT_traj = c(28, 30, 5, 15),
                     RT_states = as.integer(1:4))
  retra <- define_retra(data = data, d = d,
                        trajectories = trajectories, states = states)

  data2 <- c("28[1-1]", "30[2-2]", "5[3-3]", "15[4-4]")
  retra2 <- define_retra(data = data2, d = d,
                        trajectories = trajectories, states = states)

  expect_equal(retra, retra2)

  expect_equal(length(retra$newT$Segments), nrow(data))
  expect_equal(retra$newT$Size, nrow(data))
  expect_equal(retra$newT$Length, sum(d["28_1", "30_2"], d["30_2", "5_3"],
                                   d["5_3", "15_4"]))
  expect_equal(sum(retra$newT$Link_distance$Distance), retra$newT$Length)

})

test_that("defines circular trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  dimnames(d) <- list(paste0(trajectories, "_", states),
                      paste0(trajectories, "_", states))

  data <- data.frame(RT = rep("newT", 6),
                     RT_traj = c(28, 30, 5, 15, 28, 28),
                     RT_states = as.integer(c(1:4, 1, 4)))
  retra <- define_retra(data = data, d = d,
                        trajectories = trajectories, states = states)

  data2 <- c("28[1-1]", "30[2-2]", "5[3-3]", "15[4-4]", "28[1-1]", "28[4-4]")
  retra2 <- define_retra(data = data2, d = d,
                         trajectories = trajectories, states = states)

  expect_equal(retra, retra2)

  expect_equal(length(retra$newT$Segments), nrow(data))
  expect_equal(retra$newT$Size, nrow(unique(data[, c("RT_traj", "RT_states")])))
  expect_equal(retra$newT$Length, sum(d["28_1", "30_2"], d["30_2", "5_3"],
                                      d["5_3", "15_4"], d["15_4", "28_1"],
                                      d["28_1", "28_4"]))
  expect_equal(sum(retra$newT$Link_distance$Distance), retra$newT$Length)

})

test_that("returns errors", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  retra <- retra_edr(d, trajectories = trajectories, states = states, minSegs = 5)

  data_incomplete <- data.frame(RT = rep("A", 6),
                                RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))
  data_numericstates <- data.frame(RT = rep("A", 6),
                                   RT_traj = c(1, 1, 1, 2, 3, 3),
                                   RT_states = c(1, 2, 3, 3, 4, 5))
  data_retra <- data.frame(RT = rep("A", 6),
                     RT_traj = c(1, 1, 1, 2, 3, 3),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))

  d_df <- data.frame(d)
  short_traj <- trajectories[2:length(trajectories)]
  short_states <- states[2:length(states)]

  expect_error(define_retra(data_incomplete),
               regexp = "'data' must contain at least three columns named \"RT\", \"RT_traj\", and \"RT_states\".")
  expect_error(define_retra(data = data_numericstates),
               regexp = "The column \"RT_states\" in 'data' needs to be of class integer.")
  expect_warning(define_retra(data = data_retra, retra = retra),
                 regexp = "'data' does not contain a column named \"RT_retra\".")
  expect_error(define_retra(data = data_retra, d = d_df, trajectories = trajectories, states = states),
               regexp = "'d' must be a symmetric dissimilarity matrix or an object of class 'dist'.")
  expect_error(define_retra(data = data_retra, d = d, trajectories = short_traj, states = states),
               regexp = "The length of 'trajectories' must be equal to both dimensions in 'd'.")
  expect_error(define_retra(data = data_retra, d = d, trajectories = trajectories, states = short_states),
               regexp = "The length of 'states' must be equal to both dimensions in 'd'.")


})








