test_that("dissimilarity of an EDR to itself is zero", {
  requireNamespace("vegan", quietly = TRUE)
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  trajectories <- paste0(abun$EDR, "_", abun$traj)
  states <- abun$state

  dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates,
                                                                    sites = trajectories,
                                                                    surveys = states))
  metrics <- c("dDR", "minDist", "maxDist")

  dEDR_traj <- lapply(setNames(metrics, metrics), function(imetric){
    dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
             edr = rep(1:3, each = 30),
             metric = imetric,
             symmetrize = NULL)
  })
  dEDR_states <- lapply(setNames(metrics, metrics), function(imetric){
    dist_edr(d = as.matrix(dStates), d.type = "dStates",
             trajectories = trajectories, states = states,
             edr = rep(1:3, each = 150),
             metric = imetric,
             symmetrize = NULL)
  })

  expect_equal(dEDR_traj, dEDR_states)

  expect_equal(unique(diag(dEDR_traj[["dDR"]])), 0)
  expect_equal(unique(diag(dEDR_traj[["minDist"]])), 0)
  expect_equal(unique(diag(dEDR_traj[["maxDist"]])), 0)

})

test_that("symmetrize argument works", {
  EDR4 <- EDR_data$EDR3$abundance[traj %in% 1:15]
  EDR4$EDR <- 4
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance,
                EDR4)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates,
                                                                    sites = paste0(abun$EDR, "_", abun$traj),
                                                                    surveys = abun$state))

  dEDR_asym <- dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
                        edr = c(rep(1:3, each = 30), rep(4, 15)),
                        metric = "dDR",
                        symmetrize = NULL)
  dEDR_asym_low <- dEDR_asym[lower.tri(dEDR_asym)]
  dEDR_asym_upp <- t(dEDR_asym)[lower.tri(dEDR_asym)]

  symmetrization <- c("mean", "min", "max", "lower", "upper")
  dEDR <- lapply(setNames(symmetrization, symmetrization), function(isym){
    dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
             edr = c(rep(1:3, each = 30), rep(4, 15)),
             metric = "dDR",
             symmetrize = isym)
  })

  # test results for d.type = "dStates"
  dEDR_asym_st <- dist_edr(d = as.matrix(dStates), d.type = "dStates",
                           trajectories = paste0(abun$EDR, "_", abun$traj),
                           states = abun$state,
                           edr = c(rep(1:3, each = 150), rep(4, 15*5)),
                           metric = "dDR",
                           symmetrize = NULL)
  dEDR_st <- lapply(setNames(symmetrization, symmetrization), function(isym){
    dist_edr(d = as.matrix(dStates), d.type = "dStates",
             trajectories = paste0(abun$EDR, "_", abun$traj),
             states = abun$state,
             edr = c(rep(1:3, each = 150), rep(4, 15*5)),
             metric = "dDR",
             symmetrize = isym)
  })
  expect_equal(dEDR_asym_st, dEDR_asym)
  expect_equal(dEDR_st, dEDR)


  # Test symmetry
  expect_true(isSymmetric(dEDR[["mean"]]))
  expect_true(isSymmetric(dEDR[["min"]]))
  expect_true(isSymmetric(dEDR[["max"]]))
  expect_true(isSymmetric(dEDR[["lower"]]))
  expect_true(isSymmetric(dEDR[["upper"]]))

  # Test values
  expect_equal(dEDR[["mean"]][lower.tri(dEDR[["mean"]])],
               vapply(1:length(dEDR_asym_low), function(x){
                 mean(c(dEDR_asym_low[x], dEDR_asym_upp[x]))
               }, numeric(1)))

  expect_equal(dEDR[["min"]][lower.tri(dEDR[["min"]])],
               vapply(1:length(dEDR_asym_low), function(x){
                 min(c(dEDR_asym_low[x], dEDR_asym_upp[x]))
               }, numeric(1)))

  expect_equal(dEDR[["max"]][lower.tri(dEDR[["max"]])],
               vapply(1:length(dEDR_asym_low), function(x){
                 max(c(dEDR_asym_low[x], dEDR_asym_upp[x]))
               }, numeric(1)))

  expect_equal(dEDR[["lower"]][lower.tri(dEDR[["lower"]])],
               dEDR_asym_low)

  expect_equal(dEDR[["upper"]][lower.tri(dEDR[["upper"]])],
               dEDR_asym_upp)

})

test_that("the properties of dDR are fit", {
  EDR2_bis <- EDR_data$EDR2$abundance
  EDR2_bis$EDR <- 4
  EDR3_sbst <- EDR_data$EDR3$abundance[traj %in% 1:15]
  EDR3_sbst$EDR <- 5

  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance,
                EDR2_bis, EDR3_sbst)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates,
                                                                    sites = paste0(abun$EDR, "_", abun$traj),
                                                                    surveys = abun$state))
  dEDR <- dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
                   edr = c(rep(1:4, each = 30), rep(5, 15)),
                   metric = "dDR",
                   symmetrize = NULL)
  dEDR_st <- dist_edr(d = as.matrix(dStates), d.type = "dStates",
                      trajectories = paste0(abun$EDR, "_", abun$traj),
                      states = abun$state,
                      edr = c(rep(1:4, each = 150), rep(5, 15*5)),
                      metric = "dDR",
                      symmetrize = NULL)
  expect_equal(dEDR, dEDR_st)

  expect_equal(dEDR[2, 4], dEDR[4, 2])
  expect_equal(dEDR[2, 4], 0)
  expect_equal(dEDR[5, 3], 0)
  expect_true(dEDR[5, 3] < dEDR[3, 5])
  expect_true(dEDR[1, 5] >= dEDR[1, 3] &&
                dEDR[2, 5] >= dEDR[2, 3])
})

test_that("trajectories can be disordered", {
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance)
  abun <- abun[order(abun$sp1), ]
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  trajectories <- paste0(abun$EDR, "_", abun$traj)
  states <- abun$state

  dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates,
                                                                    sites = trajectories,
                                                                    surveys = states))

  dEDR <- dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
                   edr = unique(abun[, c("EDR", "traj")])$EDR,
                   metric = "dDR",
                   symmetrize = NULL)
  dEDR_st <- dist_edr(d = as.matrix(dStates), d.type = "dStates",
                      trajectories = paste0(abun$EDR, "_", abun$traj),
                      states = abun$state,
                      edr = abun$EDR,
                      metric = "dDR",
                      symmetrize = NULL)

  abun_or <- rbind(EDR_data$EDR1$abundance,
                   EDR_data$EDR2$abundance,
                   EDR_data$EDR3$abundance)
  dStates_or <- vegan::vegdist(abun_or[, -c(1:3)])
  trajectories_or <- paste0(abun_or$EDR, "_", abun_or$traj)
  states_or <- abun_or$state

  dTraj_or <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates_or,
                                                                       sites = trajectories_or,
                                                                       surveys = states_or))

  dEDR_or <- dist_edr(d = as.matrix(dTraj_or), d.type = "dTraj",
                      edr = unique(abun_or[, c("EDR", "traj")])$EDR,
                      metric = "dDR",
                      symmetrize = NULL)
  dEDR_st_or <- dist_edr(d = as.matrix(dStates_or), d.type = "dStates",
                         trajectories = paste0(abun_or$EDR, "_", abun_or$traj),
                         states = abun_or$state,
                         edr = abun_or$EDR,
                         metric = "dDR",
                         symmetrize = NULL)

  expect_equal(dEDR, dEDR_st)
  expect_equal(dEDR_or, dEDR_st_or)
  expect_equal(dEDR[as.character(1:3), as.character(1:3)], dEDR_st_or)

})

test_that("returns errors", {
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  trajectories <- paste0(abun$EDR, "_", abun$traj)
  states <- abun$state

  dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates,
                                                                    sites = trajectories,
                                                                    surveys = states))

  expect_error(dist_edr(d = as.data.frame(as.matrix(dTraj)), d.type = "dTraj",
                        edr = rep(1:3, each = 30)),
               regexp = "symmetric dissimilarity matrix")
  expect_error(dist_edr(d = as.matrix(dStates), d.type = "dStates",
                        trajectories = NULL, states = states,
                        edr = rep(1:3, each = 150)),
               regexp = "you must provide a value for 'trajectories'")
  expect_error(dist_edr(d = as.matrix(dStates), d.type = "dStates",
                        trajectories = trajectories, states = NULL,
                        edr = rep(1:3, each = 150)),
               regexp = "you must provide a value for 'states'")
  expect_error(dist_edr(d = as.matrix(dStates), d.type = "dStates",
                        trajectories = trajectories[1:5], states = states,
                        edr = rep(1:3, each = 150)),
               regexp = "The length of 'trajectories'")
  expect_error(dist_edr(d = as.matrix(dStates), d.type = "dStates",
                        trajectories = trajectories, states = states[1:5],
                        edr = rep(1:3, each = 150)),
               regexp = "The length of 'states'")
  expect_error(dist_edr(d = as.matrix(dTraj), d.type = "dTraj",
                        edr = rep(1:3, each = 10)),
               regexp = "'edr' needs to have a length")

})



