test_that("dDis, dBD, and dEve return a number in [0, 1]", {
  dTraj <- EDR_data$EDR1$traj_dissim
  dDis <- dDis_edr(d = dTraj, d.type = "dTraj",
                   trajectories = labels(dTraj),
                   reference = "1")

  dBD <- dBD_edr(d = dTraj, d.type = "dTraj",
                   trajectories = labels(dTraj))

  dEve <- dEve_edr(d = dTraj, d.type = "dTraj",
                   trajectories = labels(dTraj))

  expect_length(dDis, 1)
  expect_length(dBD, 1)
  expect_length(dEve, 1)

  expect_gte(dDis, 0)
  expect_gte(dBD, 0)
  expect_gte(dEve, 0)

  expect_lte(dDis, 1)
  expect_lte(dBD, 1)
  expect_lte(dEve, 1)

})

test_that("dDis is smaller then the reference trajectory belongs to the EDR
          than it it is external", {
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance[1:5, traj := 31][1:5, ])
  dStates <- vegan::vegdist(abun[, -c(1:3)], method = "bray")

  dDis_inEDR <- dDis_edr(d = dStates, d.type = "dStates",
                         trajectories = abun$traj, states = abun$state,
                         reference = "1")
  dDis_outEDR <- dDis_edr(d = dStates, d.type = "dStates",
                         trajectories = abun$traj, states = abun$state,
                         reference = "31")

  expect_lte(dDis_inEDR, dDis_outEDR)

})

test_that("dDis decreases when the weight of trajectories close to the reference increases",{
  dTraj <- as.matrix(EDR_data$EDR1$traj_dissim)

  closest <- which(dTraj[-1, 1] <= mean(dTraj[-1, 1]))
  w_close <- rep(1, 29)
  w_close[closest] <- 10

  farthest <- which(dTraj[-1, 1] > mean(dTraj[-1, 1]))
  w_far <- rep(1, 29)
  w_far[farthest] <- 10

  dDis_ew <- dDis_edr(d = dTraj, d.type = "dTraj", trajectories = labels(dTraj)[[1]],
                      reference = "1")

  dDis_close <- dDis_edr(d = dTraj, d.type = "dTraj", trajectories = labels(dTraj)[[1]],
                         reference = "1",
                         w.type = "precomputed", w.values = w_close)

  dDis_far <- dDis_edr(d = dTraj, d.type = "dTraj", trajectories = labels(dTraj)[[1]],
                       reference = "1",
                       w.type = "precomputed", w.values = w_far)

  expect_lte(dDis_close, dDis_ew)
  expect_lte(dDis_ew, dDis_far)

})

test_that("dBD is greater when trajectories belong to different EDRs", {
  abun <- rbind(EDR_data$EDR1$abundance[traj %in% 1:10],
                EDR_data$EDR2$abundance[traj %in% 11:20],
                EDR_data$EDR3$abundance[traj %in% 21:30])
  dStates <- vegan::vegdist(abun[, -c(1:3)], method = "bray")

  dBD_sameEDR <- dBD_edr(d = as.matrix(EDR_data$EDR1$state_dissim),
                         d.type = "dStates",
                         trajectories = EDR_data$EDR1$abundance$traj,
                         states = EDR_data$EDR1$abundance$state)
  dBD_diffEDR <- dBD_edr(d = as.matrix(dStates),
                         d.type = "dStates",
                         trajectories = abun$traj,
                         states = abun$state)

  expect_gte(dBD_diffEDR, dBD_sameEDR)

})

test_that("dEve is smaller when trajectories belong to different EDRs", {
  abun <- rbind(EDR_data$EDR1$abundance[traj %in% 1:15],
                EDR_data$EDR2$abundance[traj %in% 16:30])
  dStates <- vegan::vegdist(abun[, -c(1:3)], method = "bray")

  dEve_sameEDR <- dEve_edr(d = as.matrix(EDR_data$EDR1$state_dissim),
                           d.type = "dStates",
                           trajectories = EDR_data$EDR1$abundance$traj,
                           states = EDR_data$EDR1$abundance$state)
  dEve_diffEDR <- dEve_edr(d = as.matrix(dStates),
                           d.type = "dStates",
                           trajectories = abun$traj,
                           states = abun$state)

  expect_lte(dEve_diffEDR, dEve_sameEDR)

})

test_that("dEve is smaller when trajectories of the same EDR have greater weight", {
  abun <- rbind(EDR_data$EDR1$abundance[traj %in% 1:15],
                EDR_data$EDR2$abundance[traj %in% 16:30])
  dStates <- vegan::vegdist(abun[, -c(1:3)], method = "bray")

  dEve_ew <- dEve_edr(d = as.matrix(dStates),
                           d.type = "dStates",
                           trajectories = abun$traj,
                           states = abun$state)
  dEve_gw <- dEve_edr(d = as.matrix(dStates),
                      d.type = "dStates",
                      trajectories = abun$traj,
                      states = abun$state,
                      w.type = "precomputed",
                      w.values = c(rep(1, 15), rep(10, 15)))

  expect_lte(dEve_gw, dEve_ew)

})
