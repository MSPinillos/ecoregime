test_that("plots predicted trajectories from PETRA", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 40), sp2 = c(7, 7), sp3 = c(6, 7), sp4 = c(8, 8),
                                       sp5 = c(5, 4), sp6 = c(2, 3), sp7 = c(6, 5), sp8 = c(3, 4),
                                       sp9 = c(7, 6), sp10 = c(8, 7), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target", "target"),
                     states = as.integer(c(EDR_var$state, 1:2)),
                     targets = "target",
                     k = 5L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)
  expect_warning(plot.PETRA(x = petra,
                            petra.colors = "royalblue",
                            target.colors = "red",
                            traj.colors = "lightblue",
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "One predicted trajectory"),
                 regexp = "Trajectories will be displayed in an ordination space")
})

test_that("plots multiple predicted trajectories", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 10), sp2 = c(7, 13), sp3 = c(6, 15), sp4 = c(8, 15),
                                       sp5 = c(5, 6), sp6 = c(2, 2), sp7 = c(6, 8), sp8 = c(3, 4),
                                       sp9 = c(7, 9), sp10 = c(8, 12), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1, 1)),
                     targets = c("target1", "target2"),
                     k = 5L, minPts = 2L,
                     d_function = d_function,
                     d_args = d_args,
                     return_args = T)
  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "Multiple predicted trajectories"),
                 regexp = "Trajectories will be displayed in an ordination space")

})

test_that("plots states depending on uncertainty metrics", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 10), sp2 = c(7, 13), sp3 = c(6, 15), sp4 = c(8, 15),
                                       sp5 = c(5, 6), sp6 = c(2, 2), sp7 = c(6, 8), sp8 = c(3, 4),
                                       sp9 = c(7, 9), sp10 = c(8, 12), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1, 1)),
                     targets = c("target1", "target2"),
                     k = 10L, minPts = 2L,
                     d_function = d_function,
                     d_args = d_args,
                     return_args = T)
  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            uncert.metric = "N",
                            uncert.colors = c("pink", "purple", "black"),
                            uncert.range = c(2, 10),
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "State uncertainty"),
                 regexp = "Trajectories will be displayed in an ordination space")

})

test_that("plots target when it is the first/last state in the predicted trajectory", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = 30, sp2 = 7, sp3 = 6, sp4 = 8,
                                       sp5 = 5, sp6 = 2, sp7 = 6, sp8 = 3,
                                       sp9 = 7, sp10 = 8, sp11 = 2, sp12 = 3)
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, "target1"),
                      states = as.integer(c(EDR_var$state, 1)),
                      targets = "target1",
                      k = 10L, minPts = 2L, direction = -1,
                      d_function = d_function, d_args = d_args,
                      return_args = T)
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, "target1"),
                      states = as.integer(c(EDR_var$state, 1)),
                      targets = "target1",
                      k = 10L, minPts = 2L, direction = 1,
                      d_function = d_function, d_args = d_args,
                      return_args = T)
  expect_warning(plot.PETRA(x = petra1,
                            petra.colors = "royalblue",
                            target.colors = "red",
                            traj.colors = "lightblue",
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "target: last state"),
                 regexp = "Trajectories will be displayed in an ordination space")
  expect_warning(plot.PETRA(x = petra2,
                            petra.colors = "royalblue",
                            target.colors = "red",
                            traj.colors = "lightblue",
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "target: 1st state"),
                 regexp = "Trajectories will be displayed in an ordination space")

})

test_that("returns the same results with/without coord", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 10), sp2 = c(7, 13), sp3 = c(6, 15), sp4 = c(8, 15),
                                       sp5 = c(5, 6), sp6 = c(2, 2), sp7 = c(6, 8), sp8 = c(3, 4),
                                       sp9 = c(7, 9), sp10 = c(8, 12), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1, 1)),
                     targets = c("target1", "target2"),
                     k = 10L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  state_var_petra <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], petra$state_var))
  d <- dist(state_var_petra)
  coord <- smacof::mds(d, ndim = nrow(state_var) - 1)$conf

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "coord = NULL"),
                 regexp = "Trajectories will be displayed in an ordination space")

  plot.PETRA(x = petra,
             petra.colors = c("royalblue", "darkgreen"),
             target.colors = c("red", "orange"),
             traj.colors = "lightblue",
             coord = coord,
             trajectories = c(EDR_var$traj, petra$trajectories),
             states = as.integer(c(EDR_var$state, petra$states)),
             axes = c(1, 2),
             main = "coord is provided")

})

test_that("returns the same results when coord is not in order", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 10), sp2 = c(7, 13), sp3 = c(6, 15), sp4 = c(8, 15),
                                       sp5 = c(5, 6), sp6 = c(2, 2), sp7 = c(6, 8), sp8 = c(3, 4),
                                       sp9 = c(7, 9), sp10 = c(8, 12), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1, 1)),
                     targets = c("target1", "target2"),
                     k = 10L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  state_var_petra <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], petra$state_var))
  d <- dist(state_var_petra)
  coord <- smacof::mds(d, ndim = nrow(state_var) - 1)$conf
  order2 <- sample(1:nrow(state_var_petra))
  coord <- coord[order2, ]

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "In order"),
                 regexp = "Trajectories will be displayed in an ordination space")

  plot.PETRA(x = petra,
             petra.colors = c("royalblue", "darkgreen"),
             target.colors = c("red", "orange"),
             traj.colors = "lightblue",
             coord = coord,
             trajectories = c(EDR_var$traj, petra$trajectories)[order2],
             states = as.integer(c(EDR_var$state, petra$states))[order2],
             axes = c(1, 2),
             main = "Random")

})

test_that("returns the same results when coord includes more trajectories", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 10), sp2 = c(7, 13), sp3 = c(6, 15), sp4 = c(8, 15),
                                       sp5 = c(5, 6), sp6 = c(2, 2), sp7 = c(6, 8), sp8 = c(3, 4),
                                       sp9 = c(7, 9), sp10 = c(8, 12), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1, 1)),
                     targets = c("target1", "target2"),
                     k = 10L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  state_var_petra_plus <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], petra$state_var,
                                           EDR_data$EDR2$abundance[traj == 1, -c("EDR", "traj", "state")]))
  d <- dist(state_var_petra_plus)
  coord <- smacof::mds(d, ndim = nrow(state_var) - 1)$conf

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "coord = NULL"),
                 regexp = "Trajectories will be displayed in an ordination space")

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = c("red", "orange"),
                            traj.colors = "lightblue",
                            coord = coord,
                            trajectories = c(EDR_var$traj, petra$trajectories, rep(31, 5)),
                            states = as.integer(c(EDR_var$state, petra$states, 1:5)),
                            axes = c(1, 2),
                            main = "coord includes additional trajectories"),
                 regexp = "Only the trajectories used to compute 'x'")

})

test_that("returns the same traj.colors regardless target.colors", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 40), sp2 = c(7, 7), sp3 = c(6, 7), sp4 = c(8, 8),
                                       sp5 = c(5, 4), sp6 = c(2, 3), sp7 = c(6, 5), sp8 = c(3, 4),
                                       sp9 = c(7, 6), sp10 = c(8, 7), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target", "target"),
                     states = as.integer(c(EDR_var$state, 1:2)),
                     targets = "target",
                     k = 5L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = "royalblue",
                            target.colors = NULL,
                            traj.colors = c(rep("lightblue", 27), rep("orange", 3)),
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "Target = NULL; three orange trajs"))

  expect_warning(plot.PETRA(x = petra,
                            petra.colors = "royalblue",
                            target.colors = "red",
                            traj.colors = c(rep("lightblue", 27), rep("orange", 3)),
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "Target = red; three orange trajs"))

})


test_that("One traj color is used when the number of colors provided is wrong", {
  skip_on_cran()

  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 40, 20), sp2 = c(7, 7, 5), sp3 = c(6, 7, 7), sp4 = c(8, 8, 9),
                                       sp5 = c(5, 4, 5), sp6 = c(2, 3, 2), sp7 = c(6, 5, 7), sp8 = c(3, 4, 3),
                                       sp9 = c(7, 6, 8), sp10 = c(8, 7, 9), sp11 = c(2, 3, 2), sp12 = c(3, 4, 3))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target1", "target1", "target2"),
                     states = as.integer(c(EDR_var$state, 1:2, 1)),
                     targets = c("target1", "target2"),
                     k = 5L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  expect_warning(expect_warning(plot.PETRA(x = petra,
                                           petra.colors = c("royalblue", "darkgreen"),
                                           target.colors = "red",
                                           traj.colors = c("lightblue", "pink", "lightgreen"),
                                           uncert.metric = NULL,
                                           uncert.colors = NULL,
                                           uncert.range = NULL,
                                           coord = NULL,
                                           trajectories = NULL,
                                           states = NULL,
                                           axes = c(1, 2),
                                           main = "Only the first color in traj.colors"),
                                regexp = "The length of 'traj.colors'"))

  expect_warning(expect_warning(plot.PETRA(x = petra,
                                           petra.colors = c("royalblue", "darkgreen", "purple"),
                                           target.colors = "red",
                                           traj.colors = c("lightblue"),
                                           uncert.metric = NULL,
                                           uncert.colors = NULL,
                                           uncert.range = NULL,
                                           coord = NULL,
                                           trajectories = NULL,
                                           states = NULL,
                                           axes = c(1, 2),
                                           main = "Only the first color in petra.colors"),
                                regexp = "The length of 'petra.colors'"))

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

  class(state_var) <- c("CAP", "list")
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean",
                      d_function = "vegclust::vegdiststruct", d_args = list(x = state_var),
                      return_args = TRUE)

  expect_warning(plot.PETRA(x = petra1,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = "red",
                            traj.colors = c("lightblue"),
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "state_var is a list (cba)"))

  expect_warning(plot.PETRA(x = petra2,
                            petra.colors = c("royalblue", "darkgreen"),
                            target.colors = "red",
                            traj.colors = c("lightblue"),
                            uncert.metric = NULL,
                            uncert.colors = NULL,
                            uncert.range = NULL,
                            coord = NULL,
                            trajectories = NULL,
                            states = NULL,
                            axes = c(1, 2),
                            main = "state_var is a list (vegclust)"))

})

test_that("returns errors", {
  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.table::data.table(sp1 = c(30, 40), sp2 = c(7, 7), sp3 = c(6, 7), sp4 = c(8, 8),
                                       sp5 = c(5, 4), sp6 = c(2, 3), sp7 = c(6, 5), sp8 = c(3, 4),
                                       sp9 = c(7, 6), sp10 = c(8, 7), sp11 = c(2, 3), sp12 = c(3, 4))
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
  d_function = "dist"
  d_args = list(x = state_var)
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target", "target"),
                     states = as.integer(c(EDR_var$state, 1:2)),
                     targets = "target",
                     k = 5L, minPts = 2L,
                     d_function = d_function, d_args = d_args,
                     return_args = T)

  state_var_petra <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], petra$state_var))
  d <- dist(state_var_petra)
  coord <- smacof::mds(d, ndim = nrow(state_var) - 1)$conf

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = "median",
                          uncert.colors = NULL,
                          uncert.range = NULL,
                          coord = NULL,
                          trajectories = NULL,
                          states = NULL,
                          axes = c(1, 2)),
               regexp = "'arg' should be one of")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = "mean_dist",
                          uncert.colors = NULL,
                          uncert.range = NULL,
                          coord = NULL,
                          trajectories = NULL,
                          states = NULL,
                          axes = c(1, 2)),
               regexp = "Provide at least two non-NA values in 'uncert.colors'.")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = "mean_dist",
                          uncert.colors = c("black", "grey"),
                          uncert.range = c(1.1*min(petra$predicted_dist$mean_dist), max(petra$predicted_dist$mean_dist)),
                          coord = NULL,
                          trajectories = NULL,
                          states = NULL,
                          axes = c(1, 2)),
               regexp = "The minimum value in 'uncert.range'")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = "mean_dist",
                          uncert.colors = c("black", "grey"),
                          uncert.range = c(min(petra$predicted_dist$mean_dist), 0.9*max(petra$predicted_dist$mean_dist)),
                          coord = NULL,
                          trajectories = NULL,
                          states = NULL,
                          axes = c(1, 2)),
               regexp = "The maximum value in 'uncert.range'")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = NULL,
                          uncert.colors = NULL,
                          uncert.range = NULL,
                          coord = coord,
                          trajectories = NULL,
                          states = NULL,
                          axes = c(1, 2)),
               regexp = "Provide values for")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = NULL,
                          uncert.colors = NULL,
                          uncert.range = NULL,
                          coord = coord[6:nrow(coord), ],
                          trajectories = EDR_var$traj[6:nrow(coord)],
                          states = EDR_var$state[6:nrow(coord)]),
               regexp = "'coord' must contain coordinates for all states")

  expect_error(plot.PETRA(x = petra,
                          petra.colors = "royalblue",
                          target.colors = "red",
                          traj.colors = "lightblue",
                          uncert.metric = NULL,
                          uncert.colors = NULL,
                          uncert.range = NULL,
                          coord = coord,
                          trajectories = EDR_var$traj[6:nrow(coord)],
                          states = EDR_var$state[6:nrow(coord)]),
               regexp = "The length of 'trajectories'")

})
