test_that("returns an object of class 'RETRA'", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(1, 1, 1, 2, 3, 3),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))
  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  expect_s3_class(new_retra, "RETRA")
  expect_equal(attributes(new_retra$A)$names,
               c("minSegs", "Segments", "Size", "Length", "Link_distance", "Seg_density"))

  expect_type(new_retra$A$Segments, "character")
  expect_type(new_retra$A$Size, "integer")
  expect_type(new_retra$A$Length, "double")
  expect_type(new_retra$A$Link_distance$Link, "character")
  expect_type(new_retra$A$Link_distance$Distance, "double")

})

test_that("returns same results when data is defined from 'retra'", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  old_retra <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)

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

    data <- data.frame(RT_traj = RT_traj, RT_states)
    data$RT <- names(old_retra)[itraj]
    data$RT_retra <- names(old_retra)[itraj]

    return(data)
  })

  data <- data.frame(data.table::rbindlist(data.ls))

  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states,
                            retra = old_retra)

  expect_equal(
    lapply(old_retra, function(itraj){
      itraj <- itraj[c("minSegs", "Size", "Length", "Seg_density")]
    }),
    lapply(new_retra, function(itraj){
      itraj <- itraj[c("minSegs", "Size", "Length", "Seg_density")]
    })
  )

  expect_equal(
    lapply(old_retra, function(itraj){
      sum(itraj[["Link_distance"]]$Distance)
    }),
    lapply(new_retra, function(itraj){
      sum(itraj[["Link_distance"]]$Distance)
    })
  )

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

    return(data)
  })

  data <- data.frame(data.table::rbindlist(data.ls))

  new_retra <- define_retra(data = data, d = as.matrix(d),
                            trajectories = trajectories, states = states,
                            retra = old_retra)

  size_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Size - new_retra[[iRT]]$Size
  })
  length_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Length - new_retra[[iRT]]$Length
  })
  link_diff <- lapply(seq_along(old_retra), function(iRT){
    max(new_retra[[iRT]]$Link_distance$Distance) - max(old_retra[[iRT]]$Link_distance$Distance)
  })

  expect_true(all(size_diff >= 0))
  expect_true(all(length_diff >= 0))
  expect_true(all(link_diff >= 0))

})

test_that("the function works for some specific states", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(1, 1, 1, 2, 3, 3),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))

  expected_segments <- c("1[1-2]", "1[2-3]", "2[3-3]", "3[4-5]")
  expected_size <- 6
  expected_length <- sum(c(d[1, 2], d[2, 3], d[3, 8],
                           d[8, 8], d[8, 14], d[14, 15]))
  expected_links <- data.frame(
    Link = vapply(1:(length(expected_segments) - 1), function(i){
      paste0(expected_segments[i], " - ", expected_segments[i+1])
    }, character(1)),
    Distance = c(d[2, 2], d[3, 8], d[8, 14])
  )



  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states,
                            retra = NULL)

  expect_equal(new_retra$A$Segments, expected_segments)
  expect_equal(new_retra$A$Size, expected_size)
  expect_equal(new_retra$A$Link_distance$Link, expected_links$Link)

  expect_equal(new_retra$A$Link_distance, expected_links)
  expect_equal(new_retra$A$Link_distance$Distance, expected_links$Distance)


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








