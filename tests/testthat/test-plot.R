test_that("plots representative trajectories from RETRA-EDR", {
  skip_on_cran()
  # Empty test that returns Skip and Warning
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)

  plot(retra, d = d, trajectories = trajectories, states = states)

})

test_that("plots a fraction of retra", {
  skip_on_cran()
  # Empty test that returns Skip and Warning
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
  skip_on_cran()
  # Empty test that returns Skip and Warning
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
  skip_on_cran()
  # Empty test that returns Skip and Warning
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
  skip_on_cran()
  # Empty test that returns Skip and Warning
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
  skip_on_cran()
  # Empty test that returns Skip and Warning
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  data <- data.frame(RT = rep("A", 6),
                     RT_traj = c(28, 30, 5, 15, 28, 28),
                     RT_states = as.integer(c(1:4, 1, 4)))

  retra <- define_retra(data = data)

  plot(x = retra, d = d, trajectories = trajectories, states = states)

})

test_that("plot trajectories with different colors", {
  skip_on_cran()
  # Empty test that returns Skip and Warning
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)

  plot(retra, d = d, trajectories = trajectories, states = states,
       traj.colors = 1:length(unique(trajectories)))
})

test_that("plots selected RT, trajectories, and links in one color (each)", {
  skip_on_cran()
  # Empty test that returns Skip and Warning
  d <- as.matrix(EDR_data$EDR1$state_dissim)
  trajectories <- EDR_data$EDR1$abundance$traj
  states <- EDR_data$EDR1$abundance$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)

  plot(retra, d = d, trajectories = trajectories, states = states,
       select_RT = "T3", sel.color = "orange",
       traj.colors = "lightblue",
       link.color = "red", link.lty = NULL)

})

test_that("plots trajectories when states are not in order", {
  skip_on_cran()
  # Empty test that returns Skip and Warning
  order2 <- sample(1:nrow(EDR_data$EDR1$abundance), nrow(EDR_data$EDR1$abundance), replace = F)
  data <- EDR_data$EDR1$abundance[order2, ]
  d <- as.matrix(vegan::vegdist(data[, paste0("sp", 1:12)], method = "bray"))
  trajectories <- data$traj
  states <- data$state

  retra <- retra_edr(d = d, trajectories = trajectories,
                     states = states, minSegs = 5)

  plot(retra, d = d, trajectories = trajectories, states = states,
       select_RT = "T3", sel.color = "orange",
       traj.colors = "lightblue",
       link.color = "red", link.lty = NULL,
       main = "States are not in order")
})

