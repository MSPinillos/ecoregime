test_that("plots representative trajectories from RETRA-EDR", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots a fraction of retra", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 4),
                     RT_traj = rep(c(28, 30), each = 2),
                     RT_states = as.integer(c(1, 2, 2, 3)))

  retra <- define_retra(data = data)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots a segment", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 2),
                     RT_traj = rep(28, 2),
                     RT_states = as.integer(c(1, 2)))

  retra <- define_retra(data = data)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots a segment composed of states from two trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 2),
                     RT_traj = c(28, 30),
                     RT_states = as.integer(c(1, 2)))

  retra <- define_retra(data = data)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots a sequence of states from different trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 4),
                     RT_traj = c(28, 30, 5, 15),
                     RT_states = as.integer(1:4))

  retra <- define_retra(data = data)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots states leading to circular trajectories", {
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(28, 30, 5, 15, 28, 28),
                     RT_states = as.integer(c(1:4, 1, 4)))

  retra <- define_retra(data = data)

  plot(x = retra, d = d, trajectories = trajectories, states = states)

})

