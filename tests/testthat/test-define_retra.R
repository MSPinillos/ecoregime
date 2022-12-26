test_that("returns same results when data is defined from 'retra'", {
  d = EDR_data$EDR1$state_dissim
  trajectories = EDR_data$EDR1$abundance$traj
  states = EDR_data$EDR1$abundance$state
  old_retra <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)

  data.ls <- lapply(seq_along(old_retra), function(itraj){
    Segments <- old_retra[[itraj]]$Segments
    seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", Segments)), "-")
    RT_traj <- vapply(seg_components, function(iseg){
      c(iseg[1], iseg[1])
    }, character(2))
    RT_traj <- c(RT_traj)
    RT_states <- vapply(seg_components, function(iseg){
      c(as.integer(iseg[2]), as.integer(iseg[3]))
    }, integer(2))
    RT_states <- c(RT_states)

    data <- data.frame(RT_traj = RT_traj, RT_states)
    data$RT <- names(old_retra)[itraj]
    data$RT_retra <- names(old_retra)[itraj]

    return(data)
  })

  data <- data.frame(data.table::rbindlist(data.ls))

  new_retra <- define_retra(data = data, d = d,
                            trajectories = trajectories, states = states,
                            retra = old_retra)

  identical(lapply(old_retra, function(itraj){
    itraj <- itraj[c("minSegs", "Size", "Length", "Link_distance", "Seg_density")]
  }),
  lapply(new_retra, function(itraj){
    itraj <- itraj[c("minSegs", "Size", "Length", "Link_distance", "Seg_density")]
  }))

})
