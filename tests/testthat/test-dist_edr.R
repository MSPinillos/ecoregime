test_that("dissimilarity of an EDR to itself is zero", {
  requireNamespace("vegan", quietly = TRUE)
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  dTraj <- ecotraj::trajectoryDistances(d = dStates,
                                        sites = paste0(abun$EDR, "_", abun$traj),
                                        surveys = abun$state)
  metrics <- c("DDR", "minDist", "maxDist")
  dEDR <- lapply(setNames(metrics, metrics), function(imetric){
    dist_edr(dTraj = as.matrix(dTraj),
             edr = rep(1:3, each = 30),
             metric = imetric,
             symmetrize = NULL)
  })

  expect_equal(unique(diag(dEDR[["DDR"]])), 0)
  expect_equal(unique(diag(dEDR[["minDist"]])), 0)
  expect_equal(unique(diag(dEDR[["maxDist"]])), 0)

})

test_that("symmetrize argument works", {
  EDR4 <- EDR_data$EDR3$abundance[traj %in% 1:15]
  EDR4$EDR <- 4
  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance,
                EDR4)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  dTraj <- ecotraj::trajectoryDistances(d = dStates,
                                        sites = paste0(abun$EDR, "_", abun$traj),
                                        surveys = abun$state)
  dEDR_asym <- dist_edr(dTraj = as.matrix(dTraj),
                        edr = c(rep(1:3, each = 30), rep(4, 15)),
                        metric = "DDR",
                        symmetrize = NULL)
  dEDR_asym_low <- dEDR_asym[lower.tri(dEDR_asym)]
  dEDR_asym_upp <- t(dEDR_asym)[lower.tri(dEDR_asym)]

  symmetrization <- c("mean", "min", "max", "lower", "upper")
  dEDR <- lapply(setNames(symmetrization, symmetrization), function(isym){
    dist_edr(dTraj = as.matrix(dTraj),
             edr = c(rep(1:3, each = 30), rep(4, 15)),
             metric = "DDR",
             symmetrize = isym)
  })

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

test_that("the properties of DDR are fit", {
  EDR2_bis <- EDR_data$EDR2$abundance
  EDR2_bis$EDR <- 4
  EDR3_sbst <- EDR_data$EDR3$abundance[traj %in% 1:15]
  EDR3_sbst$EDR <- 5

  abun <- rbind(EDR_data$EDR1$abundance,
                EDR_data$EDR2$abundance,
                EDR_data$EDR3$abundance,
                EDR2_bis, EDR3_sbst)
  dStates <- vegan::vegdist(abun[, -c(1:3)])
  dTraj <- ecotraj::trajectoryDistances(d = dStates,
                                        sites = paste0(abun$EDR, "_", abun$traj),
                                        surveys = abun$state)
  dEDR <- dist_edr(dTraj = as.matrix(dTraj),
                   edr = c(rep(1:4, each = 30), rep(5, 15)),
                   metric = "DDR",
                   symmetrize = NULL)

  expect_equal(dEDR[2, 4], dEDR[4, 2])
  expect_equal(dEDR[2, 4], 0)
  expect_equal(dEDR[5, 3], 0)
  expect_true(dEDR[5, 3] < dEDR[3, 5])
  expect_true(dEDR[1, 5] >= dEDR[1, 3] &&
                dEDR[2, 5] >= dEDR[2, 3])
})




