test_that("works when trajectories and states are NOT in order (trajectories)", {
  skip_on_cran()
  st1 <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                    y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                    traj = c(rep(1:3, each = 3), 3),
                    states = c(rep(1:3, 3), 4))
  set.seed(23)
  st2 <- st1[sample(1:10, 10, replace = F), ]

  plot_edr(x = st1[, 1:2], trajectories = st1$traj, states = as.integer(st1$states),
           traj.colors = c("red", "green", "blue"), state.colors = NULL,
           type = "trajectories", main = 1)

  plot_edr(x = st2[, 1:2], trajectories = st2$traj, states = as.integer(st2$states),
           traj.colors = c( "blue", "red", "green"), state.colors = NULL,
           type = "trajectories", main = 2)

})

test_that("works when trajectories and states are NOT in order (states)", {
  skip_on_cran()
  st1 <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                    y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                    traj = c(rep(1:3, each = 3), 3),
                    states = c(rep(1:3, 3), 4))
  set.seed(23)
  st2 <- st1[sample(1:10, 10, replace = F), ]


  plot_edr(x = st1[, 1:2], trajectories = st1$traj, states = as.integer(st1$states),
           traj.colors = c("red", "green", "blue"),
           state.colors = c("yellow", "orange", "red", "grey", "green", "darkgreen",
                            "pink", "purple", "blue", "black"),
           type = "states", initial = TRUE, main = 1)

  plot_edr(x = st2[, 1:2], trajectories = st2$traj, states = as.integer(st2$states),
           traj.colors = c( "blue", "red", "green"),
           state.colors = c("purple", "red", "black", "yellow", "blue", "green",
                            "orange", "pink", "darkgreen", "grey"),
           type = "states", initial = TRUE, main = 2)

})

test_that("works when trajectories and states are NOT in order (gradient)", {
  skip_on_cran()
  st1 <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                    y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                    traj = c(rep(1:3, each = 3), 3),
                    states = c(rep(1:3, 3), 4))
  set.seed(23)
  st2 <- st1[sample(1:10, 10, replace = F), ]


  plot_edr(x = st1[, 1:2], trajectories = st1$traj, states = as.integer(st1$states),
           traj.colors = NULL, state.colors = viridis::viridis(3),
           type = "gradient", variable = st1$x, main = 1)

  plot_edr(x = st2[, 1:2], trajectories = st2$traj, states = as.integer(st2$states),
           traj.colors = NULL, state.colors = viridis::viridis(3),
           type = "gradient", variable = st2$x, main = 2)

})

test_that("works when x is dist", {
  skip_on_cran()
  st <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                   y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                   traj = c(rep(1:3, each = 3), 3),
                   states = c(rep(1:3, 3), 4))

  expect_warning(plot_edr(x = dist(st[, 1:2]), trajectories = st$traj, states = as.integer(st$states),
           traj.colors = c("red", "green", "blue"), state.colors = NULL,
           type = "trajectories", main = "trajectories"))

  expect_warning(plot_edr(x = dist(st[, 1:2]), trajectories = st$traj, states = as.integer(st$states),
                          traj.colors = c("red", "green", "blue"),
                          state.colors = c("yellow", "orange", "red", "grey", "green", "darkgreen",
                                           "pink", "purple", "blue", "black"),
                          type = "states", main = "states"))
  expect_warning(plot_edr(x = dist(st[, 1:2]), trajectories = st$traj, states = as.integer(st$states),
                          traj.colors = NULL, state.colors = viridis::viridis(3),
                          type = "gradient", variable = st$x, main = "gradient"))
})

test_that("uses one color for all trajectories and states", {
  skip_on_cran()
  st <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                   y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                   traj = c(rep(1:3, each = 3), 3),
                   states = c(rep(1:3, 3), 4))
  plot_edr(x = st[, 1:2], trajectories = st$traj, states = as.integer(st$states),
           traj.colors = "blue", state.colors = "red", initial = T,
           type = "states")
})

test_that("plots individual states", {
  skip_on_cran()
  st <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4, 0.5),
                   y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1, 2),
                   traj = c(rep(1:3, each = 3), 3, 4),
                   states = c(rep(1:3, 3), 4, 1))
  plot_edr(x = st[, 1:2], trajectories = st$traj, states = as.integer(st$states),
           traj.colors = "blue", state.colors = c(rep("blue", 10), "red"),
           type = "states")
})

test_that("returns errors", {
  st <- data.frame(x = c(0, 1, 2, 1, 2, 3, 4, 4, 4, 4),
                   y = c(0, 0, 0, 1, 2, 3, 4, 3, 2, 1),
                   traj = c(rep(1:3, each = 3), 3),
                   states = c(rep(1:3, 3), 4))

  expect_error(plot_edr(x = st$x, trajectories = st$traj, states = as.integer(st$states),
                        traj.colors = "blue", state.colors = "red",
                        type = "states"),
               regexp = "'x' needs to be of any of these classes: 'matrix', 'data.frame', 'dist'.")

  expect_error(plot_edr(x = st[, 1:2], trajectories = c(st$traj, 99), states = as.integer(st$states),
                        traj.colors = "blue", state.colors = "red",
                        type = "states"),
               regexp = "The length of 'trajectories' must be equal to the number of rows of 'x'.")

  expect_error(plot_edr(x = st[, 1:2], trajectories = st$traj, states = as.integer(c(st$states, 99)),
                        traj.colors = "blue", state.colors = "red",
                        type = "states"),
               regexp = "The length of 'trajectories' must be equal to the length of 'states'.")

  expect_error(plot_edr(x = st[, 1:2], trajectories = st$traj, states = st$states,
                        traj.colors = "blue", state.colors = "red",
                        type = "states"),
               regexp = "'states' needs to be of class integer.")

  expect_error(plot_edr(x = st[, 1:2], trajectories = st$traj, states = as.integer(st$states),
                        traj.colors = NULL, state.colors = viridis::viridis(3),
                        type = "gradient", variable = c(st$x, 10)))

  expect_error(plot_edr(x = st[, 1:2], trajectories = st$traj, states = as.integer(st$states),
                        traj.colors = NULL, state.colors = viridis::viridis(3),
                        type = "gradient", variable = c(st$x[-1], NA)),
               regexp = "There are missing values in 'variable'")

})
