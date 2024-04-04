test_that("returns expected results when method = 'nearest_state' and 'd' is not metric", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(vegan::vegdist(coord))
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  expected_values <- list(
    A = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                   reference = rep(names(reference), length(disturbed_trajectories)),
                   A_abs = c(min(d[9,1:4]) - min(d[8,1:4]),
                             min(d[9,5:7]) - min(d[8,5:7]),
                             min(d[12,1:4]) - min(d[11,1:4]),
                             min(d[12,5:7]) - min(d[11,5:7]),
                             min(d[15,1:4]) - min(d[14,1:4]),
                             min(d[15,5:7]) - min(d[14,5:7]),
                             min(d[18,1:4]) - min(d[17,1:4]),
                             min(d[18,5:7]) - min(d[17,5:7])),
                   A_rel = c((min(d[9,1:4]) - min(d[8,1:4])) / d[8,9],
                             (min(d[9,5:7]) - min(d[8,5:7])) / d[8,9],
                             (min(d[12,1:4]) - min(d[11,1:4])) / d[11,12],
                             (min(d[12,5:7]) - min(d[11,5:7])) / d[11,12],
                             (min(d[15,1:4]) - min(d[14,1:4])) / d[15,14],
                             (min(d[15,5:7]) - min(d[14,5:7])) / d[15,14],
                             (min(d[18,1:4]) - min(d[17,1:4])) / d[17,18],
                             (min(d[18,5:7]) - min(d[17,5:7])) / d[17,18])),

    Rc = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = 3,
                    reference = rep(names(reference), length(disturbed_trajectories)),
                    Rc_abs = c(min(d[9,1:4]) - min(d[10,1:4]),
                               min(d[9,5:7]) - min(d[10,5:7]),
                               min(d[12,1:4]) - min(d[13,1:4]),
                               min(d[12,5:7]) - min(d[13,5:7]),
                               min(d[15,1:4]) - min(d[16,1:4]),
                               min(d[15,5:7]) - min(d[16,5:7]),
                               min(d[18,1:4]) - min(d[19,1:4]),
                               min(d[18,5:7]) - min(d[19,5:7])),
                    Rc_rel = c((min(d[9,1:4]) - min(d[10,1:4])) / d[9,10],
                               (min(d[9,5:7]) - min(d[10,5:7])) / d[9,10],
                               (min(d[12,1:4]) - min(d[13,1:4])) / d[12,13],
                               (min(d[12,5:7]) - min(d[13,5:7])) / d[12,13],
                               (min(d[15,1:4]) - min(d[16,1:4])) / d[15,16],
                               (min(d[15,5:7]) - min(d[16,5:7])) / d[15,16],
                               (min(d[18,1:4]) - min(d[19,1:4])) / d[18,19],
                               (min(d[18,5:7]) - min(d[19,5:7])) / d[18,19])),

    NC = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = 3,
                    reference = rep(names(reference), length(disturbed_trajectories)),
                    NC_abs = c(min(d[10,1:4]) - min(d[8,1:4]),
                               min(d[10,5:7]) - min(d[8,5:7]),
                               min(d[13,1:4]) - min(d[11,1:4]),
                               min(d[13,5:7]) - min(d[11,5:7]),
                               min(d[16,1:4]) - min(d[14,1:4]),
                               min(d[16,5:7]) - min(d[14,5:7]),
                               min(d[19,1:4]) - min(d[17,1:4]),
                               min(d[19,5:7]) - min(d[17,5:7])),
                    NC_rel = c((min(d[10,1:4]) - min(d[8,1:4])) / d[8,10],
                               (min(d[10,5:7]) - min(d[8,5:7])) / d[8,10],
                               (min(d[13,1:4]) - min(d[11,1:4])) / d[11,13],
                               (min(d[13,5:7]) - min(d[11,5:7])) / d[11,13],
                               (min(d[16,1:4]) - min(d[14,1:4])) / d[14,16],
                               (min(d[16,5:7]) - min(d[14,5:7])) / d[14,16],
                               (min(d[19,1:4]) - min(d[17,1:4])) / d[17,19],
                               (min(d[19,5:7]) - min(d[17,5:7])) / d[17,19]))


  )

  calculated_values <- list(
    A = amplitude(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "nearest_state"),

    Rc = recovery(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "nearest_state"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "nearest_state")
  )

  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})

test_that("returns expected results when method = 'projection' and 'd' is metric", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(dist(coord))
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  stt <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                    target_states = 8:length(trajectories), reference = reference,
                                                    method = "projection"))

  expected_values <- list(

    Rt = data.frame(disturbed_trajectories = disturbed_trajectories,
                    Rt = c(1-d[8,9], 1-d[11,12], 1-d[14,15], 1-d[17,18])),

    A = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                   reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                   A_abs = c(NA,
                             NA,
                             NA,
                             stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance,
                             stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance,
                             stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance,
                             stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance,
                             stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance),

                   A_rel = c(NA,
                             NA,
                             NA,
                             (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance) / d[11,12],
                             (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance) / d[14,15],
                             (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance) / d[14,15],
                             (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance) / d[17,18],
                             (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance) / d[17,18])),

    Rc = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    Rc_abs = c(NA,
                               stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance,
                               NA,
                               stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance,
                               stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance,
                               stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance,
                               stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance,
                               stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance),

                    Rc_rel = c(NA,
                               (stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance) / d[9,10],
                               NA,
                               (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance) / d[12,13],
                               (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance) / d[15,16],
                               (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance) / d[15,16],
                               (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance) / d[18,19],
                               (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance) / d[18,19])),

    NC = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    NC_abs = c(- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance,
                               NA,
                               NA,
                               - stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance,
                               - stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance,
                               - stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance,
                               - stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance,
                               - stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance),

                    NC_rel = c((- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance) / d[8,10],
                               NA,
                               NA,
                               (- stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance) / d[11,13],
                               (- stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance) / d[14,16],
                               (- stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance) / d[14,16],
                               (- stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance) / d[17,19],
                               (- stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance) / d[17,19]))
  )

  calculated_values <- list(
    Rt = resistance(d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),

    A = amplitude(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "projection"),

    Rc = recovery(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "projection"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "projection")
  )

  expect_equal(expected_values$Rt, calculated_values$Rt)
  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})

test_that("returns expected results when method = 'projection' and 'd' is not metric", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(vegan::vegdist(coord))

  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  coordStates <- smacof::mds(d, ndim = nrow(d) - 1)$conf
  euc <- as.matrix(dist(coordStates))
  stt <- data.table::data.table(state_to_trajectory(d = euc, trajectories = trajectories, states = states,
                                                    target_states = as.integer(8:length(trajectories)), reference = reference,
                                                    method = "projection", coordStates = coordStates))


  expected_values <- list(

    Rt = data.frame(disturbed_trajectories = disturbed_trajectories,
                    Rt = c(1-d[8,9], 1-d[11,12], 1-d[14,15], 1-d[17,18])),

    A = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                   reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                   A_abs = c(NA,
                             NA,
                             NA,
                             NA,
                             stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance,
                             stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance,
                             stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance,
                             NA),

                   A_rel = c(NA,
                             NA,
                             NA,
                             NA,
                             (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance) / euc[14,15],
                             (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance) / euc[14,15],
                             (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance) / euc[17,18],
                             NA)),

    Rc = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    Rc_abs = c(NA,
                               NA,
                               NA,
                               stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance,
                               stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance,
                               stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance,
                               stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance,
                               NA),

                    Rc_rel = c(NA,
                               NA,
                               NA,
                               (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance) / euc[12,13],
                               (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance) / euc[15,16],
                               (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance) / euc[15,16],
                               (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance) / euc[18,19],
                               NA)),

    NC = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    NC_abs = c(- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance,
                               NA,
                               NA,
                               NA,
                               - stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance,
                               - stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance,
                               - stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance,
                               - stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance),

                    NC_rel = c((- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance) / euc[8,10],
                               NA,
                               NA,
                               NA,
                               (- stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance) / euc[14,16],
                               (- stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance) / euc[14,16],
                               (- stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance) / euc[17,19],
                               (- stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance) / euc[17,19]))
  )

  suppressWarnings(calculated_values <- list(
    Rt = resistance(d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),

    A = amplitude(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "projection"),

    Rc = recovery(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "projection"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "projection")
  ))

  expect_equal(expected_values$Rt, calculated_values$Rt)
  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})


test_that("returns expected results when method = 'mixed' and 'd' is metric", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(dist(coord))
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  stt <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                    target_states = 8:length(trajectories), reference = reference,
                                                    method = "mixed"))

  expected_values <- list(

    Rt = data.frame(disturbed_trajectories = disturbed_trajectories,
                    Rt = c(1-d[8,9], 1-d[11,12], 1-d[14,15], 1-d[17,18])),

    A = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                   reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                   A_abs = c(stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 8 & reference == "newT.1"]$distance,
                             stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 8 & reference == "newT.2"]$distance,
                             stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 11 & reference == "newT.1"]$distance,
                             stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance,
                             stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance,
                             stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance,
                             stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance,
                             stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance),

                   A_rel = c((stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 8 & reference == "newT.1"]$distance) / d[8,9],
                             (stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 8 & reference == "newT.2"]$distance) / d[8,9],
                             (stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 11 & reference == "newT.1"]$distance) / d[11,12],
                             (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance) / d[11,12],
                             (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance) / d[14,15],
                             (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance) / d[14,15],
                             (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance) / d[17,18],
                             (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance) / d[17,18])),

    Rc = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    Rc_abs = c(stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 10 & reference == "newT.1"]$distance,
                               stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance,
                               stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 13 & reference == "newT.1"]$distance,
                               stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance,
                               stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance,
                               stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance,
                               stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance,
                               stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance),

                    Rc_rel = c((stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 10 & reference == "newT.1"]$distance) / d[9,10],
                               (stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance) / d[9,10],
                               (stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 13 & reference == "newT.1"]$distance) / d[12,13],
                               (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance) / d[12,13],
                               (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance) / d[15,16],
                               (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance) / d[15,16],
                               (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance) / d[18,19],
                               (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance) / d[18,19])),

    NC = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    NC_abs = c(- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance,
                               - stt[target_state == 8 & reference == "newT.2"]$distance + stt[target_state == 10 & reference == "newT.2"]$distance,
                               - stt[target_state == 11 & reference == "newT.1"]$distance + stt[target_state == 13 & reference == "newT.1"]$distance,
                               - stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance,
                               - stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance,
                               - stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance,
                               - stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance,
                               - stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance),

                    NC_rel = c((- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance) / d[8,10],
                               (- stt[target_state == 8 & reference == "newT.2"]$distance + stt[target_state == 10 & reference == "newT.2"]$distance) / d[8,10],
                               (- stt[target_state == 11 & reference == "newT.1"]$distance + stt[target_state == 13 & reference == "newT.1"]$distance) / d[11,13],
                               (- stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance) / d[11,13],
                               (- stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance) / d[14,16],
                               (- stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance) / d[14,16],
                               (- stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance) / d[17,19],
                               (- stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance) / d[17,19]))
  )

  calculated_values <- list(
    Rt = resistance(d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),

    A = amplitude(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "mixed"),

    Rc = recovery(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "mixed"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "mixed")
  )

  expect_equal(expected_values$Rt, calculated_values$Rt)
  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})

test_that("returns expected results when method = 'mixed' and 'd' is not metric", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(vegan::vegdist(coord))

  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  coordStates <- smacof::mds(d, ndim = nrow(d) - 1)$conf
  euc <- as.matrix(dist(coordStates))
  stt <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                    target_states = as.integer(8:length(trajectories)), reference = reference,
                                                    method = "mixed", coordStates = coordStates))


  expected_values <- list(

    Rt = data.frame(disturbed_trajectories = disturbed_trajectories,
                    Rt = c(1-d[8,9], 1-d[11,12], 1-d[14,15], 1-d[17,18])),

    A = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                   reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                   A_abs = c(stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 8 & reference == "newT.1"]$distance,
                             stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 8 & reference == "newT.2"]$distance,
                             stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 11 & reference == "newT.1"]$distance,
                             stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance,
                             stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance,
                             stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance,
                             stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance,
                             stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance),

                   A_rel = c((stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 8 & reference == "newT.1"]$distance) / euc[8,9],
                             (stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 8 & reference == "newT.2"]$distance) / euc[8,9],
                             (stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 11 & reference == "newT.1"]$distance) / euc[11,12],
                             (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 11 & reference == "newT.2"]$distance) / euc[11,12],
                             (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 14 & reference == "newT.1"]$distance) / euc[14,15],
                             (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 14 & reference == "newT.2"]$distance) / euc[14,15],
                             (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 17 & reference == "newT.1"]$distance) / euc[17,18],
                             (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 17 & reference == "newT.2"]$distance) / euc[17,18])),

    Rc = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    Rc_abs = c(stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 10 & reference == "newT.1"]$distance,
                               stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance,
                               stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 13 & reference == "newT.1"]$distance,
                               stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance,
                               stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance,
                               stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance,
                               stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance,
                               stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance),

                    Rc_rel = c((stt[target_state == 9 & reference == "newT.1"]$distance - stt[target_state == 10 & reference == "newT.1"]$distance) / euc[9,10],
                               (stt[target_state == 9 & reference == "newT.2"]$distance - stt[target_state == 10 & reference == "newT.2"]$distance) / euc[9,10],
                               (stt[target_state == 12 & reference == "newT.1"]$distance - stt[target_state == 13 & reference == "newT.1"]$distance) / euc[12,13],
                               (stt[target_state == 12 & reference == "newT.2"]$distance - stt[target_state == 13 & reference == "newT.2"]$distance) / euc[12,13],
                               (stt[target_state == 15 & reference == "newT.1"]$distance - stt[target_state == 16 & reference == "newT.1"]$distance) / euc[15,16],
                               (stt[target_state == 15 & reference == "newT.2"]$distance - stt[target_state == 16 & reference == "newT.2"]$distance) / euc[15,16],
                               (stt[target_state == 18 & reference == "newT.1"]$distance - stt[target_state == 19 & reference == "newT.1"]$distance) / euc[18,19],
                               (stt[target_state == 18 & reference == "newT.2"]$distance - stt[target_state == 19 & reference == "newT.2"]$distance) / euc[18,19])),

    NC = data.frame(disturbed_trajectories = rep(disturbed_trajectories, each = length(reference)),
                    states = rep(3, length(reference)*length(disturbed_trajectories)),
                    reference = rep(c("newT.1", "newT.2"), length(disturbed_trajectories)),
                    NC_abs = c(- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance,
                               - stt[target_state == 8 & reference == "newT.2"]$distance + stt[target_state == 10 & reference == "newT.2"]$distance,
                               - stt[target_state == 11 & reference == "newT.1"]$distance + stt[target_state == 13 & reference == "newT.1"]$distance,
                               - stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance,
                               - stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance,
                               - stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance,
                               - stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance,
                               - stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance),

                    NC_rel = c((- stt[target_state == 8 & reference == "newT.1"]$distance + stt[target_state == 10 & reference == "newT.1"]$distance) / euc[8,10],
                               (- stt[target_state == 8 & reference == "newT.2"]$distance + stt[target_state == 10 & reference == "newT.2"]$distance) / euc[8,10],
                               (- stt[target_state == 11 & reference == "newT.1"]$distance + stt[target_state == 13 & reference == "newT.1"]$distance) / euc[11,13],
                               (- stt[target_state == 11 & reference == "newT.2"]$distance + stt[target_state == 13 & reference == "newT.2"]$distance) / euc[11,13],
                               (- stt[target_state == 14 & reference == "newT.1"]$distance + stt[target_state == 16 & reference == "newT.1"]$distance) / euc[14,16],
                               (- stt[target_state == 14 & reference == "newT.2"]$distance + stt[target_state == 16 & reference == "newT.2"]$distance) / euc[14,16],
                               (- stt[target_state == 17 & reference == "newT.1"]$distance + stt[target_state == 19 & reference == "newT.1"]$distance) / euc[17,19],
                               (- stt[target_state == 17 & reference == "newT.2"]$distance + stt[target_state == 19 & reference == "newT.2"]$distance) / euc[17,19]))
  )

  suppressWarnings(calculated_values <- list(
    Rt = resistance(d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),

    A = amplitude(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "mixed"),

    Rc = recovery(d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "mixed"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "mixed")
  ))

  expect_equal(expected_values$Rt, calculated_values$Rt)
  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})

test_that("returns same results when states are not in order", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13, 15, 17,
                             5, 9, 9, 13, 16,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4, 4, 5,
                             4, 10, 6, 3, 2,
                             8, 12, 12,
                             12, 20, 14))
  d <- as.matrix(vegan::vegdist(coord))

  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 5),
                    rep("T2", 5),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, 1:3, rep(1:5, 2), rep(1:3, 2)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  order2 <- sample(1:length(trajectories))
  coord2 <- coord[order2, ]
  d2 <- as.matrix(vegan::vegdist(coord2))
  trajectories2 <- trajectories[order2]
  states2 <- states[order2]

  expected_values <- list(
    Rt = resistance(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),
    A = amplitude(d = d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "nearest_state"),

    Rc = recovery(d = d, trajectories = trajectories, states = states,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "nearest_state"),

    NC = net_change(d = d, trajectories = trajectories, states = states,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "nearest_state")
  )
  calculated_values <- list(
    Rt = resistance(d = d2, trajectories = trajectories2, states = states2,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states),
    A = amplitude(d = d2, trajectories = trajectories2, states = states2,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  predisturbed_states = predisturbed_states,
                  reference = reference, method = "nearest_state"),

    Rc = recovery(d = d2, trajectories = trajectories2, states = states2,
                  disturbed_trajectories = disturbed_trajectories,
                  disturbed_states = disturbed_states,
                  reference = reference, method = "nearest_state"),

    NC = net_change(d = d2, trajectories = trajectories2, states = states2,
                    disturbed_trajectories = disturbed_trajectories,
                    disturbed_states = disturbed_states,
                    predisturbed_states = predisturbed_states,
                    reference = reference, method = "nearest_state")
  )

  expect_equal(expected_values$Rt, calculated_values$Rt)
  expect_equal(expected_values$A, calculated_values$A)
  expect_equal(expected_values$Rc, calculated_values$Rc)
  expect_equal(expected_values$NC, calculated_values$NC)

})

test_that("resistance returns errors", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- dist(coord)
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  expect_error(resistance(d = 1, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "'d' must be a symmetric dissimilarity")
  expect_error(resistance(d = d, trajectories = trajectories[-1], states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "The length of 'trajectories'")
  expect_error(resistance(d = d, trajectories = trajectories, states = states[-2],
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "The length of 'states'")
  expect_error(resistance(d = d, trajectories = trajectories, states = as.numeric(states),
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "'states' needs to be of class")
  expect_error(resistance(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = rep(disturbed_trajectories, each = 3),
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "Each value in 'disturbed_trajectories'")
  expect_error(resistance(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = "T5",
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states),
               regexp = "All 'disturbed_trajectories' must be included")
  expect_error(resistance(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = c(2,2,2,9),
                          predisturbed_states = predisturbed_states),
               regexp = "All 'disturbed_states' of a given trajectory must be included")
  expect_error(resistance(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = c(1, 1, 1, 4)),
               regexp = "All 'predisturbed_states' of a given trajectory must be included")
  expect_error(resistance(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = c(1, 1, 1, 2)),
               regexp = "A state cannot be included in both")
})

test_that("amplitude returns errors", {
  coord <- data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- dist(coord)
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  expect_error(amplitude(d = 1, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "'d' must be a symmetric dissimilarity")
  expect_error(amplitude(d = d, trajectories = paste(trajectories, "-", 1), states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "Avoid using '-'")
  expect_error(amplitude(d = d, trajectories = trajectories[-1], states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "The length of 'trajectories'")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states[-2],
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "The length of 'states'")
  expect_error(amplitude(d = d, trajectories = trajectories, states = as.numeric(states),
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "'states' needs to be of class")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = rep(disturbed_trajectories, each = 3),
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "Each value in 'disturbed_trajectories'")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = "T5",
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "All 'disturbed_trajectories' must be included")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = c(2,2,2,9),
                         predisturbed_states = predisturbed_states,
                         reference = reference),
               regexp = "All 'disturbed_states' of a given trajectory must be included")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = c(1, 1, 1, 4),
                         reference = reference),
               regexp = "All 'predisturbed_states' of a given trajectory must be included")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = c(1, 1, 1, 2),
                         reference = reference),
               regexp = "A state cannot be included in both")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = "R1"),
               regexp = "'reference' must be an object")
  expect_error(amplitude(d = d, trajectories = trajectories, states = states,
                         disturbed_trajectories = disturbed_trajectories,
                         disturbed_states = disturbed_states,
                         predisturbed_states = predisturbed_states,
                         reference = define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                                              c("R2[1-2]", "R2[2-3]", "R2[3-4]")))),
               regexp = "All states in 'reference'")
})

test_that("recovery returns errors", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- dist(coord)
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  expect_error(recovery(d = 1, trajectories = trajectories, states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "'d' must be a symmetric dissimilarity")
  expect_error(recovery(d = d, trajectories = paste(trajectories, "-", 1), states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "Avoid using '-'")
  expect_error(recovery(d = d, trajectories = trajectories[-1], states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "The length of 'trajectories'")
  expect_error(recovery(d = d, trajectories = trajectories, states = states[-2],
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "The length of 'states'")
  expect_error(recovery(d = d, trajectories = trajectories, states = as.numeric(states),
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "'states' needs to be of class")
  expect_error(recovery(d = d, trajectories = trajectories, states = states,
                        disturbed_trajectories = rep(disturbed_trajectories, each = 3),
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "Each value in 'disturbed_trajectories'")
  expect_error(recovery(d = d, trajectories = trajectories, states = states,
                        disturbed_trajectories = "T5",
                        disturbed_states = disturbed_states,
                        reference = reference),
               regexp = "All 'disturbed_trajectories' must be included")
  expect_error(recovery(d = d, trajectories = trajectories, states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = c(2,2,2,9),
                        reference = reference),
               regexp = "All 'disturbed_states' of a given trajectory must be included")
  expect_error(recovery(d = d, trajectories = trajectories, states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = "R1"),
               regexp = "'reference' must be an object")
  expect_error(recovery(d = d, trajectories = trajectories, states = states,
                        disturbed_trajectories = disturbed_trajectories,
                        disturbed_states = disturbed_states,
                        reference = define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                                             c("R2[1-2]", "R2[2-3]", "R2[3-4]")))),
               regexp = "All states in 'reference'")
})

test_that("net_change returns errors", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             5, 13, 13,
                             20, 9, 13,
                             5, 9, 9,
                             14, 16, 19,
                             14, 16, 19),
                       y = c(2, 5, 2, 2,
                             8, 8, 20,
                             4, 7, 4,
                             4, 10, 6,
                             8, 12, 12,
                             12, 20, 14))
  d <- dist(coord)
  trajectories <- c(rep("R1", 4),
                    rep("R2", 3),
                    rep("T1", 3),
                    rep("T2", 3),
                    rep("T3", 3),
                    rep("T4", 3))
  states <- as.integer(c(1:4, rep(1:3, 5)))
  reference <- define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                        c("R2[1-2]", "R2[2-3]")))

  disturbed_trajectories <- paste0("T", 1:4)
  disturbed_states <- rep(2, 4)
  predisturbed_states = rep(1, 4)

  expect_error(net_change(d = 1, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "'d' must be a symmetric dissimilarity")
  expect_error(net_change(d = d, trajectories = paste(trajectories, "-", 1), states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "Avoid using '-'")
  expect_error(net_change(d = d, trajectories = trajectories[-1], states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "The length of 'trajectories'")
  expect_error(net_change(d = d, trajectories = trajectories, states = states[-2],
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "The length of 'states'")
  expect_error(net_change(d = d, trajectories = trajectories, states = as.numeric(states),
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "'states' needs to be of class")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = rep(disturbed_trajectories, each = 3),
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "Each value in 'disturbed_trajectories'")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = "T5",
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "All 'disturbed_trajectories' must be included")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = c(2,2,2,9),
                          predisturbed_states = predisturbed_states,
                          reference = reference),
               regexp = "All 'disturbed_states' of a given trajectory must be included")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = c(1, 1, 1, 4),
                          reference = reference),
               regexp = "All 'predisturbed_states' of a given trajectory must be included")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = c(1, 1, 1, 2),
                          reference = reference),
               regexp = "A state cannot be included in both")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = "R1"),
               regexp = "'reference' must be an object")
  expect_error(net_change(d = d, trajectories = trajectories, states = states,
                          disturbed_trajectories = disturbed_trajectories,
                          disturbed_states = disturbed_states,
                          predisturbed_states = predisturbed_states,
                          reference = define_retra(data = list(c("R1[1-2]", "R1[2-3]", "R1[3-4]"),
                                                               c("R2[1-2]", "R2[2-3]", "R2[3-4]")))),
               regexp = "All states in 'reference'")
})


