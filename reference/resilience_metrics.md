# Metrics of trajectory deviation with respect to a reference trajectory

Set of metrics to analyze the deviation of disturbed trajectories from
an ecological dynamic regime (EDR) considering a representative
trajectory or a predicted trajectory as the reference. These metrics
include the resistance to the disturbance, amplitude, recovery, and net
change. Altogether, these indices inform on the resilience of a system.

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
  state_var = NULL,
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
  state_var = NULL,
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
  state_var = NULL,
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
  and `state_var` belongs.

- states:

  Vector of integers indicating the order of the states in `d` and
  `state_var` for each trajectory.

- disturbed_trajectories:

  Vector of the same class as `trajectories` indicating the identifier
  of the disturbed trajectories.

- disturbed_states:

  Vector of integers included in `states` indicating the first state
  after the release of the disturbance for each value in
  `disturbed_trajectories`.

- predisturbed_states:

  Vector of integers included in `states` indicating the last
  undisturbed state of each `disturbed_trajectories`. The previous
  states to `disturbed_states` are considered by default.

- reference:

  Object of class `RETRA` or `PETRA` indicating the trajectory taken as
  the reference to compute amplitude, recovery, and net_change of the
  disturbed trajectories (see Details). `PETRA` objects must include the
  arguments (`return_args = TRUE`).

- state_var:

  Object containing the state variables for each trajectory state
  (required if `reference` is of class `PETRA`).

- index:

  Method to calculate amplitude, recovery, or net change (`"absolute"`,
  `"relative"`; see Details).

- method:

  Method to calculate the distance between the `disturbed_states` or
  `predisturbed_states` and the `reference` trajectory. One of
  `"nearest_state"`, of `"projection"` (see Details).

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

*Amplitude* indicates the direction in which the system is deviated
during the disturbance in relation to the `reference` (Sánchez-Pinillos
et al., 2024; Sánchez-Pinillos et al., 2026).

If `reference` is of class `RETRA`: Positive values indicate that the
disturbance deviates the system towards the boundaries of the dynamic
regime. Negative values indicate that the disturbance deviates the
system towards the representative trajectory.

If `reference` is of class `PETRA`: Only values equal or greater than
zero are possible. Positive values indicate that the system deviates
from its predicted trajectory, whereas amplitude equal to zero indicates
that the system remains on its predicted trajectory.

Two indices can be calculated:

If `index = "absolute"`,

\\ A = d\_{dist,RT} - d\_{pre,RT} \\

If `index = "relative"`,

\\ A = \frac{d\_{dist,RT} - d\_{pre,RT}}{d\_{pre,dist}} \\

**Recovery (`recovery()`)**

*Recovery* quantifies the ability of the system to evolve towards the
`reference` following the relief of the disturbance (if positive) or
move in the direction of the boundaries of the dynamic regime (if
negative and `reference` is of class `RETRA`) or elsewhere (if negative
and `reference` is of class `PETRA`) (Sánchez-Pinillos et al., 2024;
Sánchez-Pinillos et al., 2026).

Two indices can be calculated:

If `index = "absolute"`,

\\ Rc = d\_{dist,RT} - d\_{post,RT} \\

If `index = "relative"`,

\\ Rc = \frac{d\_{dist,RT} - d\_{post,RT}}{d\_{dist,post}} \\

**Net change (`net_change()`)**

*Net change* quantifies the proximity of the system to the `reference`
relative to the pre-disturbed state after the disturbance
(Sánchez-Pinillos et al., 2024; Sánchez-Pinillos et al., 2026).

If reference is of class `RETRA`: Positive values indicate that the
system eventually evolves towards the boundaries of the dynamic regime.
Negative values indicate that the system eventually evolves towards the
`reference`.

If `reference` is of class `PETRA`: Only values equal or greater than
zero are possible. Positive values indicate that the system eventually
deviates from its predicted trajectory, whereas net change equal to zero
indicates that the system reaches its predicted trajectory some time
after the disturbance.

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
by two different methods:

- If `method = "nearest_state"`, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated as the dissimilarity between the
  pre-disturbance, disturbed, or post-disturbance states and their
  nearest state in the `reference`.

- If `method = "projection"`, \\d\_{pre,RT}\\, \\d\_{dist,RT}\\, and
  \\d\_{post,RT}\\ are calculated as the dissimilarity between the
  pre-disturbance, disturbed, or post-disturbance states and their
  projection onto the `reference` or using the nearest state of the
  `reference` if the dissimilarity is smaller.

## References

Sánchez-Pinillos, M., Leduc, A., Ameztegui, A., Kneeshaw, D., Lloret,
F., & Coll, L. (2019). Resistance, resilience or change:
Post-disturbance dynamics of boreal forests after insect outbreaks.
*Ecosystems* 22, 1886-1901 https://doi.org/10.1007/s10021-019-00378-6

Sánchez-Pinillos, M., Dakos, V., & Kéfi, S. (2024). Ecological dynamic
regimes: A key concept for assessing ecological resilience. *Biological
Conservation* 289, 110409 https://doi.org/10.1016/j.biocon.2023.110409

Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
Forecasting ecological trajectories from ecological dynamic regimes to
improve resilience analysis. *Methods in Ecology and Evolution*.
<doi:10.1111/2041-210x.70372>

## See also

[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
to identify representative trajectories in an ecological dynamic regime.

[`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)
to generate an object of class`RETRA`.

[`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md)
to predict ecological trajectories in an ecological dynamic regime.

[`MPD()`](https://mspinillos.github.io/ecoregime/reference/MPD.md) for
estimating the prediction accuracy of `petra-edr()` outputs.

[`state_to_trajectory()`](https://mspinillos.github.io/ecoregime/reference/state_to_trajectory.md)
to calculate the position of a state with respect to a trajectory.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE)) {
  # Calculate the resistance of disturbed systems

  # Abundance matrix including disturbed and undisturbed trajectories
  abundance <- rbind(EDR_data$EDR3$abundance,
                     EDR_data$EDR3_disturbed$abundance, fill = TRUE)

  # State dissimilarities (Bray-Curtis) for disturbed and undisturbed trajectories
  d <- vegan::vegdist(abundance[, paste0("sp", 1:12)], method = "bray")

  # Resistance
  Rt <- resistance(d = d, trajectories = abundance$traj, states = abundance$state,
                   disturbed_trajectories = unique(abundance[!is.na(disturbed_states)]$traj),
                   disturbed_states = abundance[disturbed_states == 1]$state)

  # Calculate the amplitude, recovery, and net change of disturbed systems taking
  # a representative trajectory as the reference. Analogous metrics can be
  # calculated using a PETRA object as the reference.

  # Identify the representative trajectories of the EDR from undisturbed trajectories
  RT <- retra_edr(d = EDR_data$EDR3$state_dissim,
                  trajectories = EDR_data$EDR3$abundance$traj,
                  states = as.integer(EDR_data$EDR3$abundance$state),
                  minSegs = 5)

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
