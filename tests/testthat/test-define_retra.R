test_that("returns an object of class 'RETRA'", {
  d = as.matrix(EDR_data$EDR2$state_dissim)
  trajectories = EDR_data$EDR2$abundance$traj
  states = EDR_data$EDR2$abundance$state
  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(1, 1, 1, 2, 3, 3),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))
  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  expect_s3_class(new_retra, "RETRA")
})

test_that("returns same results when data is defined from 'retra'", {
  d = EDR_data$EDR1$state_dissim
  trajectories = EDR_data$EDR1$abundance$traj
  states = EDR_data$EDR1$abundance$state
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
  d = EDR_data$EDR2$state_dissim
  trajectories = EDR_data$EDR2$abundance$traj
  states = EDR_data$EDR2$abundance$state
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

  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states,
                            retra = old_retra)

  size_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Size - new_retra[[iRT]]$Size
  })
  length_diff <- lapply(seq_along(old_retra), function(iRT){
    old_retra[[iRT]]$Length - new_retra[[iRT]]$Length
  })
  link_diff <- lapply(seq_along(old_retra), function(iRT){
    new_retra[[iRT]]$Length - max(old_retra[[iRT]]$Link_distance$Distance)
  })

  expect_true(all(size_diff > 0))
  expect_true(all(length_diff > 0))
  expect_true(all(link_diff > 0))

})

test_that("the function works for some specific states", {
  d = as.matrix(EDR_data$EDR2$state_dissim)
  trajectories = EDR_data$EDR2$abundance$traj
  states = EDR_data$EDR2$abundance$state
  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(1, 1, 1, 2, 3, 3),
                     RT_states = as.integer(c(1, 2, 3, 3, 4, 5)))
  expected_segments <- c("1[1-2]", "1[2-3]", "2[3-3]", "3[4-5]")
  expected_size <- 6
  expected_length <- sum(c(d[1, 2], d[2, 3], d[3, 8],
                           d[8, 8], d[8, 14], d[14, 15]))
  expected_links <- c(0, d[3, 8], d[8, 14])

  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states)

  expect_equal(new_retra$A$Segments, expected_segments)
  expect_equal(new_retra$A$Size, expected_size)
  expect_equal(new_retra$A$Link_distance$Distance, expected_links)

})
