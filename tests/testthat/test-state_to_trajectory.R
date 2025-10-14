test_that("returns expected results when 'd' is metric and reference is not 'RETRA'", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- dist(coord)
  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- c("T1", "T2")

  expected_values <- list(
    NS = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(reference, length(target_states)),
                    distance = c(3, sqrt(45), 3, 4, 5, 5),
                    relative_position = c(0, 0, 5/20, 0, 10/20, 8/20)),
    P = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(reference, length(target_states)),
                    distance = c(3, sqrt(45), 3, 0, 3, 5),
                    relative_position = c(0, 0, 5/20, 4/20, 14/20, 8/20))
  )

  calculated_values <- lapply(setNames(c("nearest_state", "projection"), c("NS", "P")), function(imethod){
    state_to_trajectory(d = d, trajectories = trajectories,
                        states = states, target_states = target_states,
                        reference = reference, method = imethod)
  })

  expect_equal(expected_values$NS, calculated_values$NS)
  expect_equal(expected_values$P, calculated_values$P)

})

test_that("returns same results when states are not in order", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- dist(coord)
  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- c("T1", "T2")

  set.seed(123)
  order2 <- sample(1:length(trajectories))
  coord2 <- coord[order2, ]
  d2 <- as.matrix(dist(coord2))
  trajectories2 <- trajectories[order2]
  states2 <- states[order2]
  target_states2 <- which(trajectories2 %in% LETTERS[1:3])


  expected_values <- lapply(setNames(c("nearest_state", "projection"), c("NS", "P")), function(imethod){
    df <- state_to_trajectory(d = d, trajectories = trajectories,
                              states = states, target_states = target_states,
                              reference = reference, method = imethod)
    df$target_state <- trajectories[df$target_state]
    return(df)
  })

  calculated_values <- lapply(setNames(c("nearest_state", "projection"), c("NS", "P")), function(imethod){
    df <- state_to_trajectory(d = d2, trajectories = trajectories2,
                              states = states2, target_states = target_states2,
                              reference = reference, method = imethod)
    df$target_state <- trajectories2[df$target_state]
    df <- df[order(df$target_state),]
    row.names(df) <- 1:nrow(df)
    return(df)
  })

  expect_equal(expected_values$NS, calculated_values$NS)
  expect_equal(expected_values$P, calculated_values$P)

})


test_that("returns expected results when 'd' is metric and reference is 'RETRA'", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- dist(coord)
  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- define_retra(data = list(c("T1[1-2]", "T1[2-3]", "T1[3-4]"),
                                        c("T2[1-2]", "T2[2-3]")))

  expected_values <- list(
    NS = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(names(reference), length(target_states)),
                    distance = c(3, sqrt(45), 3, 4, 5, 5),
                    relative_position = c(0, 0, 5/20, 0, 10/20, 8/20)),
    P = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(names(reference), length(target_states)),
                    distance = c(3, sqrt(45), 3, 0, 3, 5),
                    relative_position = c(0, 0, 5/20, 4/20, 14/20, 8/20))
  )

  calculated_values <- lapply(setNames(c("nearest_state", "projection"), c("NS", "P")), function(imethod){
    state_to_trajectory(d = d, trajectories = trajectories,
                        states = states, target_states = target_states,
                        reference = reference, method = imethod)
  })

  expect_equal(expected_values$NS, calculated_values$NS)
  expect_equal(expected_values$P, calculated_values$P)

})

test_that("returns expected results when 'd' is not metric and coordStates is provided", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- as.matrix(vegan::vegdist(coord))
  coordStates <- round(cmdscale(d), 3)

  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- define_retra(data = list(c("T1[1-2]", "T1[2-3]", "T1[3-4]"),
                                        c("T2[1-2]", "T2[2-3]")))

  euc <- as.matrix(dist(coordStates))


  expected_values <- list(
    NS = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(names(reference), length(target_states)),
                    distance = c(d[5,1], d[5,8], d[6,2], d[6,9], d[7,3], d[7,9]),
                    relative_position = c(0, 0,
                                          d[1,2]/sum(d[1,2], d[2,3], d[3,4]), d[8,9]/sum(d[8,9], d[9,10]),
                                          sum(d[1,2], d[2,3])/sum(d[1,2], d[2,3], d[3,4]), d[8,9]/sum(d[8,9], d[9,10]))),
    P = data.frame(target_state = rep(target_states, each = length(reference)),
                    reference = rep(names(reference), length(target_states)),
                    distance = c(euc[5,1], euc[5,8], euc[6,2],
                                 sqrt(euc[6,8]^2 - ((euc[6,8]^2 - euc[6,9]^2 + euc[8,9]^2)/(2*euc[8,9]))^2),
                                 sqrt(euc[7,2]^2 - ((euc[7,2]^2 - euc[7,3]^2 + euc[2,3]^2)/(2*euc[2,3]))^2),
                                 euc[7,9]),
                    relative_position = c(0, 0, euc[1,2]/sum(c(euc[1,2], euc[2,3], euc[3,4])),
                                          ((euc[6,8]^2 - euc[6,9]^2 + euc[8,9]^2)/(2*euc[8,9])) / sum(c(euc[8,9], euc[9,10])),
                                          sum(c(euc[1,2], (euc[7,2]^2 - euc[7,3]^2 + euc[2,3]^2)/(2*euc[2,3]))) / sum(c(euc[1,2], euc[2,3], euc[3,4])),
                                          euc[8,9]/sum(c(euc[8,9], euc[9,10]))))
  )

  calculated_values <- lapply(setNames(c("nearest_state", "projection"), c("NS", "P")), function(imethod){
    state_to_trajectory(d = d, trajectories = trajectories,
                        states = states, target_states = target_states,
                        reference = reference, method = imethod, coordStates = coordStates)
  })

  expect_equal(expected_values$NS, calculated_values$NS)
  expect_equal(expected_values$P, calculated_values$P)

})

test_that("returns expected results when 'd' is not metric and coordStates is NULL", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- as.matrix(vegan::vegdist(coord))
  coordStates <- smacof::mds(d, ndim = nrow(d) - 1)$conf

  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- c("T1", "T2")
  euc <- as.matrix(dist(coordStates))

  expected_values <- data.frame(target_state = rep(target_states, each = length(reference)),
                                reference = rep(reference, length(target_states)),
                                distance = c(euc[5,1],
                                             euc[5,8], euc[6,2],
                                             sqrt(euc[6,8]^2 - ((euc[6,8]^2 - euc[6,9]^2 + euc[8,9]^2)/(2*euc[8,9]))^2),
                                             sqrt(euc[7,3]^2 - ((euc[7,3]^2 - euc[7,4]^2 + euc[3,4]^2)/(2*euc[3,4]))^2),
                                             euc[7,9]),
                                relative_position = c(0,
                                                      0, euc[1,2]/sum(c(euc[1,2], euc[2,3], euc[3,4])),
                                                      ((euc[6,8]^2 - euc[6,9]^2 + euc[8,9]^2)/(2*euc[8,9])) / sum(c(euc[8,9], euc[9,10])),
                                                      sum(c(euc[1,2], euc[2,3], (euc[7,3]^2 - euc[7,4]^2 + euc[3,4]^2)/(2*euc[3,4]))) / sum(c(euc[1,2], euc[2,3], euc[3,4])),
                                                      euc[8,9]/sum(c(euc[8,9], euc[9,10]))))

  suppressWarnings(calculated_values <- state_to_trajectory(d = d, trajectories = trajectories,
                        states = states, target_states = target_states,
                        reference = reference, method = "projection", coordStates = NULL))

  expect_equal(expected_values, calculated_values)

})

test_that("returns same or different results depending on the number of target_states", {
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


  stt <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                    target_states = as.integer(c(9, 17)), reference = reference,
                                                    method = "projection", coordStates = NULL))
  expect_warning(stt2 <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                                    target_states = as.integer(c(8, 17)), reference = reference,
                                                                    method = "projection", coordStates = NULL)))
  stt3 <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                     target_states = as.integer(17), reference = reference,
                                                     method = "projection", coordStates = NULL))
  stt4 <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                     target_states = as.integer(c(8, 17)), reference = reference,
                                                     method = "nearest_state", coordStates = NULL))
  stt5 <- data.table::data.table(state_to_trajectory(d = d, trajectories = trajectories, states = states,
                                                     target_states = as.integer(17), reference = reference,
                                                     method = "nearest_state", coordStates = NULL))

  expect_false(isTRUE(all.equal(stt2[target_state == 17], stt3)))
  expect_equal(stt[target_state == 17], stt3)
  expect_equal(stt4[target_state == 17], stt5)

})

test_that("returns errors", {
  coord <-  data.frame(x = c(5, 9, 13, 23,
                             2, 9, 17,
                             5, 13, 13),
                       y = c(2, 5, 2, 2,
                             2, 8, 5,
                             8, 8, 20))
  d <- dist(coord)
  trajectories <- c("T1", "T1", "T1", "T1",
                    "A", "B", "C",
                    "T2", "T2", "T2")
  states <- as.integer(c(1:4, 1, 1, 1, 1:3))
  target_states <- which(trajectories %in% LETTERS[1:3])
  reference <- c("T1", "T2")
  retra <- define_retra(data = list(c("T1[1-2]", "T1[2-3]", "T1[3-4]"),
                                    c("T2[1-2]", "T2[2-3]")))

  expect_error(state_to_trajectory(d = list(d), trajectories = trajectories,
                                   states = states, target_states = target_states,
                                   reference = reference, method = "nearest_state"),
               regexp = "'d' must be a symmetric dissimilarity")
  expect_error(state_to_trajectory(d = d, trajectories = c(trajectories, 999),
                                   states = states, target_states = target_states,
                                   reference = reference, method = "nearest_state"),
               regexp = "The length of 'trajectories'")
  expect_error(state_to_trajectory(d = d, trajectories = trajectories,
                                   states = states[-1], target_states = target_states,
                                   reference = reference, method = "nearest_state"),
               regexp = "The length of 'states'")
  expect_error(state_to_trajectory(d = d, trajectories = trajectories,
                                   states = states, target_states = c("A", "B", "C"),
                                   reference = reference, method = "nearest_state"),
               regexp = "'target_states' needs to be")
  expect_error(state_to_trajectory(d = d, trajectories = paste0(trajectories),
                                   states = states, target_states = 12,
                                   reference = reference, method = "nearest_state"),
               regexp = "'target_states' needs to be")
  expect_error(state_to_trajectory(d = d, trajectories = trajectories,
                                   states = states, target_states = target_states,
                                   reference = "T5", method = "nearest_state"),
               regexp = "'reference' must be included")
  expect_error(state_to_trajectory(d = d, trajectories = trajectories,
                                   states = states, target_states = target_states,
                                   reference = define_retra(data = list(c("T1[1-2]", "T1[2-3]", "T1[3-4]", "T1[4-5]"),
                                                                        c("T2[1-2]", "T2[2-3]", "T2[3-4]"))),
                                   method = "nearest_state"),
               regexp = "All states in 'reference' must be included")
})

