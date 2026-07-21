test_that("returns an object of class 'PETRA'", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 30]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, target_var$traj),
                     states = c(EDR_var$state, target_var$state),
                     targets = unique(target_var$traj), k = 5L, minPts = 2L,
                     d_function = "dist", d_args = list(x = state_var), return_args = T)


  expect_s3_class(petra, "PETRA")
  expect_equal(attributes(petra)$names,
               c("arguments", "k_dist", "predicted_dist", "state_var", "trajectories", "states"))

  expect_s3_class(petra$k_dist, "data.frame")
  expect_equal(class(state_var), class(petra$state_var))
  expect_equal(names(petra$arguments), c("trajectories", "states", "targets",
                                         "k", "minPts", "eps",
                                         "d_function", "d_args", "args_state_var", "d",
                                         "method", "w"))

})

test_that("returns the similar results using one point or a subtrajectory passing by the point", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- lapply(c("medoid", "mean"), function(imethod){
    petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              method = imethod, w_function = NULL, w = NULL)
  })

  #----

  target_var <- EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra2 <- lapply(c("medoid", "mean"), function(imethod){
    petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              method = imethod, w_function = NULL, w = NULL)
  })

  for (i in 1:length(petra1)) {
    expect_equal(petra1[[i]]$state_var[petra1[[i]]$states < 1, ],
                 petra2[[i]]$state_var[petra2[[i]]$states < 1, ])
  }

})

test_that("returns same results for a target regardless of the number of targets", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- lapply(c("medoid", "mean"), function(imethod){
    petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              method = imethod)
  })

  #----

  target_var <- rbind(EDR_data$EDR1$abundance[traj == 29 & state == 2],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra2 <- lapply(c("medoid", "mean"), function(imethod){
    petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              method = imethod)
  })

  for (i in 1:length(petra1)) {
    expect_equal(data.table::as.data.table(petra1[[i]]$state_var),
                 data.table::as.data.table(petra2[[i]]$state_var[petra2[[i]]$trajectories %in% unique(petra1[[i]]$trajectories), ]))
  }

})

test_that("method returns different results", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "medoid",
                      d_function = "dist", d_args = list(x = state_var))
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean",
                      d_function = "dist", d_args = list(x = state_var))
  petra3 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "linear",
                      d_function = "dist", d_args = list(x = state_var))

  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra2$state_var[which(petra2$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra3$state_var[which(petra3$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra3$state_var[which(petra3$trajectories == 30), ],
                                petra2$state_var[which(petra2$trajectories == 30), ])))
})

test_that("w_function returns different results", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "linear",
                      d_function = "dist", d_args = list(x = state_var))
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "power", alpha = 2,
                      d_function = "dist", d_args = list(x = state_var))
  petra3 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "exponential", alpha = 2,
                      d_function = "dist", d_args = list(x = state_var))
  petra4 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "Gaussian", alpha = 2,
                      d_function = "dist", d_args = list(x = state_var))
  petra5 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "hyperbolic", alpha = 2,
                      d_function = "dist", d_args = list(x = state_var))
  petra6 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean", w_function = "spherical",
                      d_function = "dist", d_args = list(x = state_var))

  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra2$state_var[which(petra2$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra3$state_var[which(petra3$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra4$state_var[which(petra4$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra5$state_var[which(petra5$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra1$state_var[which(petra1$trajectories == 30), ],
                                petra6$state_var[which(petra6$trajectories == 30), ])))

  expect_false(isTRUE(all.equal(petra2$state_var[which(petra2$trajectories == 30), ],
                                petra3$state_var[which(petra3$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra2$state_var[which(petra2$trajectories == 30), ],
                                petra4$state_var[which(petra4$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra2$state_var[which(petra2$trajectories == 30), ],
                                petra5$state_var[which(petra5$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra2$state_var[which(petra2$trajectories == 30), ],
                                petra6$state_var[which(petra6$trajectories == 30), ])))

  expect_false(isTRUE(all.equal(petra3$state_var[which(petra3$trajectories == 30), ],
                                petra4$state_var[which(petra4$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra3$state_var[which(petra3$trajectories == 30), ],
                                petra5$state_var[which(petra5$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra3$state_var[which(petra3$trajectories == 30), ],
                                petra6$state_var[which(petra6$trajectories == 30), ])))

  expect_false(isTRUE(all.equal(petra4$state_var[which(petra4$trajectories == 30), ],
                                petra5$state_var[which(petra5$trajectories == 30), ])))
  expect_false(isTRUE(all.equal(petra4$state_var[which(petra4$trajectories == 30), ],
                                petra6$state_var[which(petra6$trajectories == 30), ])))

  expect_false(isTRUE(all.equal(petra5$state_var[which(petra5$trajectories == 30), ],
                                petra6$state_var[which(petra6$trajectories == 30), ])))

})

test_that("works when w is provided", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  w <- calculate_w(d = dist(state_var),
                   trajectories = c(EDR_var$traj, target_var$traj),
                   states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                   targets = unique(target_var$traj),
                   k = 5L,
                   w_function = "linear", alpha = 1)

  petra1 <- petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              w = w)

  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w_function = "linear")

  petra3 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w = rep(0, length(w)))

  petra4 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var))

  expect_equal(petra1, petra2)
  expect_equal(petra3, petra4)

  #----

  target_var <- rbind(EDR_data$EDR1$abundance[traj == 29 & state == 2],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  w <- calculate_w(d = dist(state_var),
                   trajectories = c(EDR_var$traj, target_var$traj),
                   states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                   targets = unique(target_var$traj),
                   k = 5L,
                   w_function = "linear", alpha = 1)

  petra5 <- petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              w = w)

  petra6 <- petra_edr(state_var = state_var,
              trajectories = c(EDR_var$traj, target_var$traj),
              states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
              targets = unique(target_var$traj), k = 5L, minPts = 2L,
              d_function = "dist", d_args = list(x = state_var),
              w_function = "linear")

  w[[1]] <- rep(0, length(w[[1]]))
  petra7 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w = w)

  petra8 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = as.integer(c(EDR_var$state, 1:length(target_var$traj)-1)),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var),
                      w_function = c(NA, "linear"))

  expect_equal(petra5, petra6)
  expect_equal(petra7, petra8)

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
                      d_function = "cba::sdists", d_args = list(x = state_var))

  class(state_var) <- c("CAP", "list")
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      method = "mean",
                      d_function = "vegclust::vegdiststruct", d_args = list(x = state_var))

  expect_s3_class(petra1, "PETRA")
  expect_s3_class(petra2, "PETRA")

})

test_that("returns same results when d is not null", {
  requireNamespace("vegan", quietly = TRUE)
  EDR_var <- EDR_data$EDR1$abundance[traj < 30]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  d <- vegan::vegdist(state_var, method = "bray")

  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "vegan::vegdist", d_args = list(x = state_var, method = "bray"))

  petra2 <- petra_edr(state_var = state_var, d = d,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj), k = 5L, minPts = 2L,
                      d_function = "vegan::vegdist", d_args = list(x = state_var, method = "bray"))

  expect_equal(petra1, petra2)

})

test_that("k, eps, and minPts work for one target", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 30]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  petra1 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj),
                      k = 2L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var))
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj),
                      k = 30L, minPts = 2L,
                      d_function = "dist", d_args = list(x = state_var))

  petra3 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj),
                      k = 30L, minPts = 10L,
                      d_function = "dist", d_args = list(x = state_var))

  petra4 <- petra_edr(state_var = state_var,
                      trajectories = c(EDR_var$traj, target_var$traj),
                      states = c(EDR_var$state, target_var$state),
                      targets = unique(target_var$traj),
                      k = 30L, minPts = 2L, eps = 10,
                      d_function = "dist", d_args = list(x = state_var))

  expect_lte(length(petra1$ID_petra), length(petra2$ID_petra))
  expect_lte(length(petra3$ID_petra), length(petra2$ID_petra))
  expect_lte(length(petra4$ID_petra), length(petra2$ID_petra))

})

test_that("k, eps, and minPts work for multiple targets and values", {
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
                      method = "mean", w_function = "linear", alpha = NA, w = NULL)

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
                      method = "mean", w_function = NULL, alpha = NA, w = NULL)

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

  expect_equal(data.table::data.table(petra1$state_var),
               data.table::data.table(petra4$state_var[petra4$trajectories == unique(petra1$trajectories), ]))
  expect_equal(data.table::data.table(petra2$state_var),
               data.table::data.table(petra4$state_var[petra4$trajectories == unique(petra2$trajectories), ]))
  expect_equal(data.table::data.table(petra3$state_var),
               data.table::data.table(petra4$state_var[petra4$trajectories == unique(petra3$trajectories), ]))
})

test_that("returns the same results regardless of the order of the states", {
  edr_data <- EDR_data$EDR1$abundance
  cols <- names(edr_data)[-c(1)]
  targets <- edr_data[c(1, 6, 11), ][, (cols) := lapply(.SD, "+", 2), .SDcols = cols]
  targets$traj <- 31:33

  abundance <- rbind(edr_data, targets)
  d <- dist(abundance[, paste0('sp', 1:12)])

  state_var <- data.frame(abundance[, paste0("sp", 1:12)])
  petra <- petra_edr(state_var = state_var,
                     trajectories = abundance$traj,
                     states = as.integer(abundance$state), targets = 31:33,
                     k = 150L, minPts = 2L, eps = quantile(d)[3],
                     d_function = "dist", d_args = list(x = state_var),
                     w_function = "exponential", alpha = 5, return_args = T)

  abundance2 <- abundance[sample(1:nrow(abundance)), ]
  d2 <- dist(abundance2[, paste0('sp', 1:12)])

  state_var2 <- data.frame(abundance2[, paste0("sp", 1:12)])
  petra2 <- petra_edr(state_var = state_var,
                      trajectories = abundance$traj,
                      states = as.integer(abundance$state), targets = 31:33,
                      k = 150L, minPts = 2L, eps = quantile(d)[3],
                      d_function = "dist", d_args = list(x = state_var),
                      w_function = "exponential", alpha = 5, return_args = T)

  expect_equal(petra$state_var, petra2$state_var)

})

test_that("there are no neighboring states", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 29]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))

  expect_warning(petra1 <- petra_edr(state_var = state_var,
                                     trajectories = c(EDR_var$traj, target_var$traj),
                                     states = as.integer(c(EDR_var$state, target_var$state)),
                                     targets = unique(target_var$traj),
                                     k = 5L, minPts = 2L, eps = 0.1,
                                     d_function = "dist", d_args = list(x = state_var)))

  target_var <- EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  expect_warning(petra2 <- petra_edr(state_var = state_var,
                                     trajectories = c(EDR_var$traj, target_var$traj),
                                     states = as.integer(c(EDR_var$state, target_var$state)),
                                     targets = unique(target_var$traj),
                                     k = 5L, minPts = 2L, eps = 0.1,
                                     d_function = "dist", d_args = list(x = state_var)))

  expect_equal(petra1$k_dist$k_dist, NA)
  expect_equal(petra2$k_dist$k_dist, c(NA, NA))

})

test_that("returns errors", {
  requireNamespace("vegan", quietly = TRUE)
  EDR_var <- EDR_data$EDR1$abundance[traj < 28]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  d <- dist(state_var)

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 2L, minPts = 5L,
                         d_function = "dist", d_args = list(x = state_var)),
               "'minPts' cannot be greater than 'k'.")

  expect_error(petra_edr(state_var = c(1, 3, 5),
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var)),
               "'state_var' needs to be of any of these classes: 'matrix', 'data.frame', 'list'")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = EDR_var$traj,
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var)),
               "The length of 'trajectories' must be equal to the length of 'states'.")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = as.character(c(EDR_var$state, target_var$state)),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var)),
               "'states' needs to be of class integer.")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = EDR_var$traj,
                         states = EDR_var$state,
                         targets = EDR_var$traj[1], k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var)),
               "The number of elements in 'state_var' must be equal to the length of")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = "A", k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var)),
               "'target' needs to be included in 'trajectories.")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = EDR_var$traj[1], k = 5L, minPts = 2L,
                         d_function = "vegan::vegdist", d_args = list(state_var = state_var)))

  # expect_error(petra_edr(state_var = state_var,
  #                        trajectories = c(EDR_var$traj, target_var$traj),
  #                        states = c(EDR_var$state, target_var$state),
  #                        targets = EDR_var$traj[1], k = 5L, minPts = 2L,
  #                        d_function = "vegan::vegdist", d_args = list(abun = state_var)),
  #              "'d_function' must return an object")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean", w = 1))

  expect_error(petra_edr(state_var = state_var[-1, ],
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5, minPts = 2,
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         d = list(d),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean", direction = 0))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         d = as.matrix(d)[-1, -1],
                         method = "mean"))

  #--------------

  target_var <- rbind(EDR_data$EDR1$abundance[traj == 28 & state %in% 2],
                      EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 3:4])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  d <- dist(state_var)

  w <- calculate_w(d = d, trajectories = c(EDR_var$traj, target_var$traj),
                   states = c(EDR_var$state, target_var$state),
                   targets = unique(target_var$traj), k = c(5L, 2L, 2L),
                   w_function = "linear")

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = c(5L, 2L), minPts = 2L,
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = c(5L, 2L, 2L), minPts = c(4L, 1L),
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         eps = c(0.1, 0.5),
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         w = w[1:2],
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         w = lapply(w, as.character),
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))

  expect_error(petra_edr(state_var = state_var,
                         trajectories = c(EDR_var$traj, target_var$traj),
                         states = c(EDR_var$state, target_var$state),
                         targets = unique(target_var$traj), k = 5L, minPts = 2L,
                         w_function = c("linear", NA),
                         d_function = "dist", d_args = list(x = state_var),
                         method = "mean"))



})

#----------------------------------

test_that("calculate_w returns a numeric vector of the length of trajectories when
          there is one target of one state", {
            EDR_var <- EDR_data$EDR1$abundance[traj < 27]
            target_var <- rbind(EDR_data$EDR1$abundance[traj == 27 & state == 2])
            state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                          target_var[, -c("EDR", "traj", "state")]))
            trajectories = c(EDR_var$traj, target_var$traj)
            states = as.integer(c(EDR_var$state, target_var$state))
            targets = unique(target_var$traj)

            d <- dist(state_var)

            w <- calculate_w(d = d, trajectories = trajectories, k = length(trajectories),
                             states = states, targets = targets,
                             w_function = "linear", alpha = 1)

            expect_true(inherits(w, "numeric"))
            expect_equal(length(w), length(trajectories))

          })

test_that("calculate_w returns a numeric vector of the length of trajectories when
          there is one multiple-state target", {
            EDR_var <- EDR_data$EDR1$abundance[traj < 27]
            target_var <- rbind(EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
            state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                          target_var[, -c("EDR", "traj", "state")]))
            trajectories = c(EDR_var$traj, target_var$traj)
            states = as.integer(c(EDR_var$state, target_var$state))
            targets = unique(target_var$traj)

            d <- dist(state_var)

            w <- calculate_w(d = d, trajectories = trajectories,
                             k = length(trajectories),
                             states = states, targets = targets,
                             w_function = "linear", alpha = 1)

            expect_true(inherits(w, "numeric"))
            expect_equal(length(w), length(trajectories))

          })


test_that("calculate_w returns a list of the length of targets when
          there are one-state and multiple-state targets", {
            EDR_var <- EDR_data$EDR1$abundance[traj < 27]
            target_var <- rbind(EDR_data$EDR1$abundance[traj == 27 & state == 2],
                                EDR_data$EDR1$abundance[traj == 28 & state == 2],
                                EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3],
                                EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
            state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                          target_var[, -c("EDR", "traj", "state")]))
            trajectories = c(EDR_var$traj, target_var$traj)
            states = as.integer(c(EDR_var$state, target_var$state))
            targets = unique(target_var$traj)

            d <- dist(state_var)

            w <- calculate_w(d = d, trajectories = trajectories,
                             k = length(trajectories),
                             states = states, targets = targets,
                             w_function = "power", alpha = 2)

            expect_true(inherits(w, "list"))
            expect_equal(length(w), length(targets))
            expect_equal(unique(vapply(w, length, numeric(1))),
                         length(trajectories))
          })

test_that("calculate_w assigns w=0 to one-state targets and each target itself", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 27]
  target_var <- rbind(EDR_data$EDR1$abundance[traj == 27 & state == 2],
                      EDR_data$EDR1$abundance[traj == 28 & state == 2],
                      EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  trajectories = c(EDR_var$traj, target_var$traj)
  states = as.integer(c(EDR_var$state, target_var$state))
  targets = unique(target_var$traj)

  d <- dist(state_var)

  w <- calculate_w(d = d, trajectories = trajectories,
                   states = states, targets = targets,
                   k = length(trajectories),
                   w_function = "exponential", alpha = 3)

  expect_equal(unique(as.numeric(vapply(seq_along(w), function(itarget){
    w[[itarget]][trajectories %in% 27:28]
  }, numeric(2)))), 0)
  expect_equal(unique(unlist(lapply(seq_along(w), function(itarget){
    w[[itarget]][trajectories == targets[itarget]]
  }))), 0)

})

test_that("calculate_w works for multiple targets and values", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 28]
  target_var1 <- EDR_data$EDR1$abundance[traj == 28 & state %in% 1]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var1[, -c("EDR", "traj", "state")]))
  trajectories1 <- c(EDR_var$traj, target_var1$traj)

  w1 <- calculate_w(d = dist(state_var),
                    trajectories = trajectories1,
                    states = as.integer(c(EDR_var$state, target_var1$state)),
                    targets = unique(target_var1$traj),
                    k = 5L, eps = 7,
                    w_function = "linear", alpha = NA)

  #----

  target_var2 <- EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var2[, -c("EDR", "traj", "state")]))
  trajectories2 <- c(EDR_var$traj, target_var2$traj)

  w2 <- c(rep(1, length(EDR_var$traj)),
          rep(0, length(target_var2$traj)))

  #----

  target_var3 <- EDR_data$EDR1$abundance[traj == 30 & state %in% 3:4]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var3[, -c("EDR", "traj", "state")]))
  trajectories3 <- c(EDR_var$traj, target_var3$traj)

  w3 <- calculate_w(d = dist(state_var),
                    trajectories = trajectories3,
                    states = as.integer(c(EDR_var$state, target_var3$state)),
                    targets = unique(target_var3$traj),
                    k = 4L, eps = NULL,
                    w_function = "exponential", alpha = 3)

  #----

  target_var4 <- rbind(EDR_data$EDR1$abundance[traj == 28 & state == 1],
                       EDR_data$EDR1$abundance[traj == 29 & state %in% 2:3],
                       EDR_data$EDR1$abundance[traj == 30 & state %in% 3:4])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var4[, -c("EDR", "traj", "state")]))
  trajectories4 <- c(EDR_var$traj, target_var4$traj)

  w4 <- calculate_w(d = dist(state_var),
                    trajectories = trajectories4,
                    states = as.integer(c(EDR_var$state, target_var4$state)),
                    targets = unique(target_var4$traj),
                    k = c(5L, 140L, 4L), eps = c(7, 10, NA),
                    w_function = c("linear", NA, "exponential"), alpha = c(NA, NA, 3))


  expect_equal(w1, w4[[1]][trajectories4 %in% trajectories1])
  expect_equal(w2, w4[[2]][trajectories4 %in% trajectories2])
  expect_equal(w3, w4[[3]][trajectories4 %in% trajectories3])

})

test_that("calculate_w returns errors", {
  EDR_var <- EDR_data$EDR1$abundance[traj < 28]
  target_var <- EDR_data$EDR1$abundance[traj == 30 & state == 2]
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  d <- dist(state_var)


  expect_error(calculate_w(d = list(d),
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = as.matrix(d)[-1, -1],
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = as.matrix(d),
                           trajectories = c(EDR_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.numeric(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = "A",
                           k = 5L, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5.0, eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = "5",
                           w_function = "linear", alpha = 1))

  #---------------

  target_var <- rbind(EDR_data$EDR1$abundance[traj == 28 & state == 2],
                      EDR_data$EDR1$abundance[traj == 29 & state %in% 3:4],
                      EDR_data$EDR1$abundance[traj == 30 & state %in% 2:3])
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("EDR", "traj", "state")]))
  d <- dist(state_var)

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = c(5L, 2L), eps = NULL,
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = c(7, NA),
                           w_function = "linear", alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = c("linear", "exponential"), alpha = 1))

  expect_error(calculate_w(d = d,
                           trajectories = c(EDR_var$traj, target_var$traj),
                           states = as.integer(c(EDR_var$state, target_var$state)),
                           targets = unique(target_var$traj),
                           k = 5L, eps = NULL,
                           w_function = c(NA, "exponential", "exponential"), alpha = c(1, 2)))

})

