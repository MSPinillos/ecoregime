test_that("returns a numeric vector of the same length than the number of targets", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 10]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var), return_args = T)

  mpd1 <- MPD(petra1)
  expect_equal(length(mpd1), length(unique(target_var$traj)))

  #----------

  EDR_var <- EDR_data$EDR1$abundance[traj < 10]
  target_var <- rbind(EDR_data$EDR1$abundance[traj == 28 & state == 2],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var), return_args = T)
  mpd2 <- MPD(petra2)

  expect_equal(length(mpd2), length(unique(target_var$traj)))

})

test_that("MPD is larger for targets out of the EDR", {
  state_var <- data.frame(rbind(EDR_data$EDR1$abundance[traj < 10, -c("EDR", "traj", "state")],
                                EDR_data$EDR1$abundance[traj == 10 & state == 1, -c("EDR", "traj", "state")],
                                data.frame(sp1 = 0, sp2 = 0, sp3 = 0, sp4 = 0, sp5 = 0, sp6 = 0,
                                           sp7 = 0, sp8 = 0, sp9 = 0, sp10 = 0, sp11 = 100, sp12 = 0)))
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_data$EDR1$abundance[traj < 10]$traj, 10, "target"),
                     states = as.integer(c(EDR_data$EDR1$abundance[traj < 10]$state, 1, 1)),
                     targets = c(10, "target"), d_function = "dist", d_args = list(x = state_var),
                     k = 3L, minPts = 2L, return_args = T)

  mpd <- MPD(petra)

  expect_gt(mpd[2], mpd[1])

})

test_that("the error is reduced when w != NULL", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 10]
  target_var <- EDR_data$EDR1$abundance[traj == 10 & state == 1]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w_function = NULL,
                      return_args = T)

  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w_function = "exponential", alpha = 100,
                      return_args = T)

  mpd1 <- MPD(x = petra1)
  mpd2 <- MPD(x = petra2)

  expect_gte(mpd1, mpd2)

})

test_that("works when state_var is a list", {
  skip_on_cran()
  skip_if_not_installed("cba")
  skip_if_not_installed("vegclust")

  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj %in% 29:30 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  state_var <- lapply(1:nrow(state_var), function(i) {
    state_var[i, ]
  })

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean",
                      d_function = "cba::sdists", d_args = list(x = state_var),
                      return_args = TRUE)
  expect_true(is.numeric(MPD(petra1)))

})

test_that("works when there are multiple conditions for different targets", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 28]
  target_var <- EDR_data$EDR1$abundance[traj == 28 & state %in% 1]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, target_var$state)),
                      targets = unique(target_var$traj),
                      k = 5L, minPts = 2L, eps = 7,
                      d_function = "dist", d_args = list(x = state_var),
                      method = "mean", w_function = "linear", alpha = NA, w = NULL,
                      return_args = T)
  mpd1 <- MPD(petra1)

  #----

  target_var <- EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, target_var$state)),
                      targets = unique(target_var$traj),
                      k = 140L, minPts = 5L, eps = 10,
                      d_function = "dist", d_args = list(x = state_var),
                      method = "mean", w_function = NULL, alpha = NA, w = NULL,
                      return_args = T)

  mpd2 <- MPD(petra2)

  #----

  target_var <- EDR_data$EDR1$abundance[traj == 30 & state %in% 3:4]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra3 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, target_var$state)),
                      targets = unique(target_var$traj),
                      k = 4L, minPts = 4L, eps = NULL,
                      d_function = "dist", d_args = list(x = state_var),
                      method = "mean", w_function = "exponential", alpha = 3, w = NULL,
                      return_args = T)
  mpd3 <- MPD(petra3)

  #----

  target_var <- rbind(EDR_data$EDR1$abundance[traj == 28 & state == 1],
                      EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 3:4])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra4 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, target_var$state)),
                      targets = unique(target_var$traj),
                      k = as.integer(c(5, 140, 4)), minPts = as.integer(c(2, 5, 4)),
                      eps = c(7, 10, NA),
                      d_function = "dist", d_args = list(x = state_var),
                      method = "mean", w_function = c("linear", NA, "exponential"),
                      alpha = c(NA, NA, 3), w = NULL,
                      return_args = T)

  mpd4 <- MPD(petra4)


  expect_equal(c(mpd1, mpd2, mpd3), mpd4)

})


test_that("returns errors", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 10]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      return_args = F)

  expect_error(MPD(state_var),
               "'x' must be an object of class 'PETRA'.")
  expect_error(MPD(petra1),
               "'x' must contain 'arguments'.")
})


