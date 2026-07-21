test_that("summarizes representative trajectories from retra_edr and define_retra", {
  d <- EDR_data$EDR1$state_dissim
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state
  retra <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)
  sum_retra <- summary(retra)

  data_A <- data.frame(RT = rep("A", 4), # name of the new representative trajectory
                       RT_traj = c(3, 3, 4, 4), # identifier of the original trajectories
                       RT_states = c(1:2, 4:5)) # states in the original trajectories
  data_B <- data.frame(RT = rep("B", 5), # name of the new representative trajectory
                       RT_traj = c(5, 5, 6, 7, 7), # identifier of the original trajectories
                       RT_states = c(1, 2, 4, 4, 5)) # states in the original trajectories
  df <- rbind(data_A, data_B)
  df$RT_states <- as.integer(df$RT_states)
  defined <- define_retra(data = df)
  sum_defined <- summary(defined)

  expect_equal(nrow(sum_retra), length(retra))
  expect_type(sum_retra$Size, "double")
  expect_type(sum_retra$Length, "double")
  expect_type(sum_retra$Avg_link, "double")
  expect_type(sum_retra$Sum_link, "double")
  expect_type(sum_retra$Avg_density, "double")
  expect_type(sum_retra$Max_density, "double")
  expect_type(sum_retra$Avg_depth, "double")
  expect_type(sum_retra$Max_depth, "double")

  expect_equal(nrow(sum_defined), length(defined))
  expect_type(sum_defined$Size, "double")
  expect_true(all(is.na(sum_defined$Length)))
  expect_true(all(is.na(sum_defined$Avg_link)))
  expect_true(all(is.na(sum_defined$Sum_link)))
  expect_true(all(is.na(sum_defined$Avg_density)))
  expect_true(all(is.na(sum_defined$Max_density)))
  expect_true(all(is.na(sum_defined$Avg_depth)))
  expect_true(all(is.na(sum_defined$Max_depth)))
})
