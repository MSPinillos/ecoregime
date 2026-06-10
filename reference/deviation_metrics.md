# Metrics of trajectory deviation with respect to a reference trajectory

Set of metrics to analyze the deviation of disturbed trajectories from
an ecological dynamic regime (EDR) considering a representative
trajectory as the reference. These metrics include the resistance to the
disturbance, amplitude, recovery, and net change.

## Usage

``` r
resistance(
  d,
  trajectories,
  states,
  disturbed_trajectories,
  disturbed_states,
  predisturbed_states = disturbed_states - 1
)

amplitude(
  d,
  trajectories,
  states,
  disturbed_trajectories,
  disturbed_states,
  predisturbed_states = disturbed_states - 1,
  reference,
  index = c("absolute", "relative"),
  method = "nearest_state"
)

recovery(
  d,
  trajectories,
  states,
  disturbed_trajectories,
  disturbed_states,
  reference,
  index = c("absolute", "relative"),
  method = "nearest_state"
)

net_change(
  d,
  trajectories,
  states,
  disturbed_trajectories,
  disturbed_states,
  predisturbed_states = disturbed_states - 1,
  reference,
  index = c("absolute", "relative"),
  method = "nearest_state"
)
```

## Arguments

- d:

  Either a symmetric matrix or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states.

- trajectories:

  Vector indicating the trajectory or site to which each state in `d`
  belongs.

- states:

  Vector of integers indicating the order of the states in `d` for each
  trajectory.

- disturbed_trajectories:

  Vector of the same class as `trajectories` indicating the identifier
  of the disturbed trajectories.

- disturbed_states:

  Vector of integers included in `states`indicating the first state
  after the release of the disturbance for each value in
  `disturbed_trajectories`.

- predisturbed_states:

  Vector of integers included in `states` indicating the last
  undisturbed state of each `disturbed_trajectories`. The previous
  states to `disturbed_states` are considered by default.

- reference:

  Object of class `RETRA` indicating the representative trajectory taken
  as the reference to compute the amplitude, recovery, and net_change of
  the disturbed trajectories (see Details).

- index:

  Method to calculate amplitude, recovery, or net change (`"absolute"`,
  `"relative"`; see Details).

- method:

  Method to calculate the distance between the `disturbed_states` or
  `predisturbed_states` and the `reference` trajectory. One of
  `"nearest_state"`, `"projection"` or `"mixed"` (see Details).

## Value

- `resistance()` returns a data frame of two columns indicating the
  resistance value (`Rt`) for each `disturbed_trajectory`.

- `amplitude()` returns a data frame of three columns indicating the
  amplitude value (`A_abs`; `A_rel`) for each `disturbed_trajectory` and
  `reference`. If `index = c("absolute", "relative")`, both values are
  included in a data frame of four columns.

- `recovery()` returns a data frame of four columns indicating the
  recovery value (`Rc_abs`; `Rc_rel`) for each `disturbed_trajectory`,
  post-disturbance state (`state`) and `reference`. If
  `index = c("absolute", "relative")`, both values are included in a
  data frame of five columns.

- `net_change` returns a data frame of four columns indicating the net
  change value (`NC_abs`; `NC_rel`) for each `disturbed_trajectory`,
  post-disturbance state (`state`), and `reference`. If
  `index = c("absolute", "relative")`, both values are included in a
  data frame of five columns.

## Details

**Resistance (`resistance()`)**

*Resistance* captures the immediate impact of the disturbance as a
function of the changes in the state variables (Sánchez-Pinillos et al.,
2019).

\\ Rt = 1 - d\_{pre,dist} \\

**Amplitude (`amplitude()`)**

*Amplitude* indicates the direction in which the system is displaced
during the disturbance in relation to the `reference` (Sánchez-Pinillos
et al., 2024). Positive values indicate that the disturbance displaces
the system towards the boundaries of the dynamic regime. Negative values
indicate that the disturbance displaces the system towards the
representative trajectory.

Two indices can be calculated:

If `index = "absolute"`,

\\ A = d\_{dist,RT} - d\_{pre,RT} \\

If `index = "relative"`,

\\ A = \frac{d\_{dist,RT} - d\_{pre,RT}}{d\_{pre,dist}} \\

**Recovery (`recovery()`)**

*Recovery* quantifies the ability of the system to evolve towards the
`reference` following the relief of the disturbance (if positive) or
move in the direction of the boundaries of the dynamic regime (if
negative) (Sánchez-Pinillos et al., 2024).

Two indices can be calculated:

If `index = "absolute"`,

\\ Rc = d\_{dist,RT} - d\_{post,RT} \\

If `index = "relative"`,

\\ Rc = \frac{d\_{dist,RT} - d\_{post,RT}}{d\_{dist,post}} \\

**Net change (`net_change()`)**

*Net change* quantifies the proximity of the system to the `reference`
relative to the pre-disturbed state (Sánchez-Pinillos et al., 2024).
Positive values indicate that the system eventually evolves towards the
boundaries of the dynamic regime. Negative values indicate that the
system eventually evolves towards the `reference`.

Two indices can be calculated:

If `index = "absolute"`,

\\ NC = d\_{post,RT} - d\_{pre,RT} \\

If `index = "relative"`,

\\ NC = \frac{d\_{post,RT} - d\_{pre,RT}}{d\_{pre,post}} \\

In all cases:

- \\d\_{pre,RT}\\ is the dissimilarity between the `predisturbed_states`
  and the `reference`.

- \\d\_{dist,RT}\\ is the dissimilarity between the `disturbed_states`
  and the `reference`.

- \\d\_{post,RT}\\ is the dissimilarity between the states after
  `disturbed_states` and the `reference`.

- \\d\_{pre,dist}\\ is the dissimilarity contained in `d` between the
  `predisturbed_states` and the `disturbed_states`.

- \\d\_{dist,post}\\ is the dissimilarity contained in `d` between the
  `disturbed_states` and the post-disturbed states.

- \\d\_{pre,post}\\ is the dissimilarity contained in `d` between the
  `predisturbed_states` and the post-disturbed states.

\\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and \\d\_{post,RT}\\ are calculated
using the function
[`state_to_trajectory()`](https://mspinillos.github.io/ecoregime/reference/state_to_trajectory.md)
by three different methods:

- If `method = "nearest_state"`, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated as the dissimilarity between the
  pre-disturbance, disturbed, or post-disturbance states and their
  nearest state in the `reference`.

- If `method = "projection"`, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated as the dissimilarity between the
  pre-disturbance, disturbed, or post-disturbance states and their
  projection onto the `reference`.

- If `method = "mixed"`, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated in the same way than
  `method = "projection"` whenever the pre-disturbance, disturbed and
  post-disturbance states can be projected onto any segment of the
  `reference`. Otherwise, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated using the nearest state of the
  `reference`.

## References

Sánchez-Pinillos, M., Leduc, A., Ameztegui, A., Kneeshaw, D., Lloret,
F., & Coll, L. (2019). Resistance, resilience or change:
Post-disturbance dynamics of boreal forests after insect outbreaks.
*Ecosystems* 22, 1886-1901 https://doi.org/10.1007/s10021-019-00378-6

Sánchez-Pinillos, M., Dakos, V., & Kéfi, S. (2024). Ecological dynamic
regimes: A key concept for assessing ecological resilience. *Biological
Conservation* 289, 110409 https://doi.org/10.1016/j.biocon.2023.110409

## See also

[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
to identify representative trajectories in an ecological dynamic regime.

[`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)
to generate an object of class`RETRA`.

[`state_to_trajectory()`](https://mspinillos.github.io/ecoregime/reference/state_to_trajectory.md)
to calculate the position of a state with respect to a trajectory.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE)) {
  # Identify the representative trajectories of the EDR from undisturbed trajectories
  RT <- retra_edr(d = EDR_data$EDR3$state_dissim,
                  trajectories = EDR_data$EDR3$abundance$traj,
                  states = as.integer(EDR_data$EDR3$abundance$state),
                  minSegs = 5)

  # Abundance matrix including disturbed and undisturbed trajectories
  abundance <- rbind(EDR_data$EDR3$abundance,
                     EDR_data$EDR3_disturbed$abundance, fill = TRUE)

  # State dissimilarities (Bray-Curtis) for disturbed and undisturbed trajectories
  d <- vegan::vegdist(abundance[, paste0("sp", 1:12)], method = "bray")

  # Resistance
  Rt <- resistance(d = d, trajectories = abundance$traj, states = abundance$state,
                   disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
                   disturbed_states = abundance[disturbed_states == 1]$state)

  # Amplitude
  A <- amplitude(d = d, trajectories = abundance$traj, states = abundance$state,
                 disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
                 disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)

  # Recovery
  Rc <- recovery(d = d, trajectories = abundance$traj, states = abundance$state,
                 disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
                 disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)

  # Net change
  NC <- net_change(d = d, trajectories = abundance$traj, states = abundance$state,
                   disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
                   disturbed_states = abundance[disturbed_states == 1]$state, reference = RT)
}
```
