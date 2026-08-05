# Predicted Ecological Trajectories in Ecological Dynamic Regimes (PETRA-EDR)

Predicts ecological trajectories from a target state or a sequence of
states based on their k-nearest states in the trajectories forming a
reference ecological dynamic regime (EDR) and the previous and following
states in their respective trajectories.

## Usage

``` r
petra_edr(
  state_var,
  trajectories,
  states,
  targets,
  d_function,
  d_args,
  d = NULL,
  k,
  eps = NULL,
  minPts = k,
  w_function = NULL,
  alpha = 1,
  w = NULL,
  method = "mean",
  direction = 2,
  return_args = FALSE
)
```

## Arguments

- state_var:

  Object containing the state variables for each trajectory state in
  both the reference dynamic regime and the `targets`. The class must
  match the one required in `d_function`. It can be a list if that is
  required in `d_function`.

- trajectories:

  Vector indicating the trajectory or site to which each state in
  `state_var` belongs.

- states:

  Vector of integers indicating the order of the states in `state_var`
  for each trajectory.

- targets:

  Vector indicating the trajectory or site of the target or targets for
  which longer trajectories must be predicted.

- d_function:

  Either a function or a non-empty character string naming the function
  to be called to compute state dissimilarities (see
  [`do.call`](https://rdrr.io/r/base/do.call.html)). The output of
  `d_function` must be an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html).

- d_args:

  A list of arguments to the `d_function` call (see
  [`do.call`](https://rdrr.io/r/base/do.call.html)).

- d:

  Either a symmetric matrix or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states in `state_var`. The
  elements need to follow the order indicated by `trajectories` and
  `states`. If `NULL`, `d` is calculated by calling `d_function`.

- k:

  Vector of integers indicating the number of near states to the target
  in the trajectories composing the reference ecological dynamic regime.

- eps:

  Numeric vector indicating the dissimilarity threshold for each target
  beyond which trajectories in the dynamic regime are not considered.

- minPts:

  Vector of integers including the minimum number of states required to
  predict preceding and following states for each target.

- w_function:

  String indicating the weighting function to estimate the state
  variables of the predicted states for each target: `"linear"`,
  `"power"`, `"exponential"`, `"Gaussian"`, `"hyperbolic"`,
  `"spherical"`.

- alpha:

  Numeric value of the shape parameter used by `w_function` (required
  for `"power"`, `"exponential"`, and `"hyperbolic"`).

- w:

  Numeric vector of the same length as `trajectories` giving the weights
  to estimate the state variables in the predicted trajectory. If
  `targets` is a vector of two or more elements, list of the same length
  than `targets` containing a numeric vector with trajectory weights for
  each target.

- method:

  String naming the method used to compute the states forming the
  predicted trajectory: `"mean"` (default) or `"medoid"`.

- direction:

  String or integer indicating whether the prediction is done backward
  (`"backward"` or `-1`), forward (`"forward"` or `1`), or in both
  directions (`"both"` or `2`, default).

- return_args:

  Logic value indicating whether a list containing the arguments
  provided in the function call must be returned.

## Value

The function `petra_edr()` returns an object of class `PETRA`, which is
a list with elements containing the following information:

- `arguments`:

  List of arguments as specified when the function is called (only if
  `return_args = TRUE`).

- `k_dist`:

  Data frame of five columns including `target`: identifier of the
  target provided in the argument `targets`; `target_state`: indices of
  the states in the extremes of the `targets`; `k_trajectories`: vector
  indicating the trajectory or site to which each of the k-nearest
  states belongs; `k_states`: vector indicating the order of the
  k-nearest states in their trajectory (`k_trajectories`); `k_dist`:
  dissimilarity between the target and each of the k-nearest states.

- `predicted_dist`:

  Data frame of six columns including, for each predicted state
  (`predicted_state`) of the target (`target`) the number of states in
  the EDR used in the averaging process to calculate the predicted
  states (`N`; \\N \ge MinPts\\) and the mean (`mean_dist`), standard
  deviation (`sd_dist`), minimum (`min_dist`), and maximum (`max_dist`)
  dissimilarities to the predicted states.

- `state_var`:

  Updated `state_var` object containing the state variables of the
  states forming the predicted trajectory/ies.

- `trajectories`:

  Vector indicating the target to which each state in the predicted
  `state_var` object belongs.

- `states`:

  Vector of integers indicating the order of the states in the predicted
  `state_var` object for each target. In order to keep the original
  state value of the targets, there can be values lower or equal to
  zero.

## Details

The PETRA-EDR algorithm (`petra_edr()`) is based on the search for the
`k` nearest states of the trajectories forming a reference EDR to the
target and a moving window towards both directions of their
corresponding trajectories.

First, PETRA-EDR looks for the target's `k` nearest states in the
trajectories forming the reference EDR. Alternatively, PETRA-EDR may
consider all states in a pre-defined radius `eps` by assigning `k` the
total number of states in the trajectories of the reference EDR.

Then, consecutive states forming the predicted trajectory are defined
based on the previous and following states of the `k` nearest states in
their respective trajectories.

Predicted states are defined by averaging the values of the state
variables of at least `minPts` states in the trajectories to which the
k-nearest states belong (`method = "mean"`) or using the state variables
of the medoid of those states (`method = "medoid"`).

To estimate a weighted mean of the state variables, it is necessary to
indicate the weight (`w`) that must be used in each trajectory of the
EDR. Different weights can be used for each target by providing a list
of numeric vectors with the weights assigned to each trajectory state.
Alternatively, one can specify the weighting function to be applied
depending on the dissimilarity between the target and the states in the
trajectories forming the EDR:

**`"linear"`**

\$\$w(d_i) = 1 - \frac{d_i}{d\_{max}}\$\$

**`"power"`**

\$\$w(d_i) = 1 - {(\frac{d_i}{d\_{max}})^{\alpha}}\$\$

**`"exponential"`**

\$\$w(d_i) = e^{\frac{- \alpha d_i}{d\_{max}}}\$\$

**`"Gaussian"`**

\$\$w(d_i) = e^{-(\frac{d_i}{d\_{max}})^2}\$\$

**`"hyperbolic"`**

\$\$w(d_i) = \frac{(1 + \frac{d_i}{d\_{max}})^{- \alpha} - 2^{-
\alpha}}{1 - 2^{- \alpha}}\$\$

**`"spherical"`**

\$\$w(d_i) = 1 - 1.5 \frac{d_i}{d\_{max}} +
0.5(\frac{d_i}{d\_{max}})^3\$\$

where \\d_i\\ is the dissimilarity between the target and the trajectory
\\i\\, \\d\_{max}\\ is the maximum dissimilarity between the target and
all trajectories, and \\\alpha\\ is a shape parameter.

## References

Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
Forecasting ecological trajectories from ecological dynamic regimes to
improve resilience analysis. *Methods in Ecology and Evolution*.
<doi:10.1111/2041-210x.70372>

## See also

[`MPD()`](https://mspinillos.github.io/ecoregime/reference/MPD.md) for
estimating the prediction accuracy of `petra-edr()` outputs.

[`plot.PETRA()`](https://mspinillos.github.io/ecoregime/reference/plot.PETRA.md)
for plotting predicted trajectories in an ordination space representing
the state space of the EDR.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE)) {
  # Example 1 -----------------------------------------------------------------
  # Compute the predicted trajectory of a target composed of one state

  # State variables for the states in the trajectories forming the EDR
  EDR_var <- EDR_data$EDR1$abundance

  # State variables of the target
  target_var <- data.table::data.table(sp1 = 30, sp2 = 10, sp3 = 13, sp4 = 12,
                                       sp5 = 5, sp6 = 2, sp7 = 6, sp8 = 3,
                                       sp9 = 7, sp10 = 8, sp11 = 2, sp12 = 3)

  # Define state_var including the state variables for the states in the EDR and
  # the target
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))

  # Define the function used to calculate state dissimilarities and its arguments
  # For example, the Canberra dissimilarity.
  d_function = "vegan::vegdist"
  d_args = list(x = state_var, method = "canberra")

  # Compute PETRA-EDR
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target"),
                     states = as.integer(c(EDR_var$state, 1)),
                     targets = "target",
                     k = 5L,
                     minPts = 2L,
                     d_function = d_function,
                     d_args = d_args)

  # Example 2 -----------------------------------------------------------------
  # Compute the predicted trajectory of two targets using different parameter
  # values

  # State variables for the states in the trajectories forming the EDR
  EDR_var <- EDR_data$EDR1$abundance

  # State variables of the target states
  target_var <- data.table::data.table(
    traj = c("target1", "target1", "target2"), state = c(1, 2, 1),
    sp1 = c(30, 31, 3), sp2 = c(10, 9, 70), sp3 = c(13, 8, 3), sp4 = c(12, 9, 4),
    sp5 = c(5, 4, 4), sp6 = c(2, 3, 3), sp7 = c(6, 7, 2), sp8 = c(3, 5, 2),
    sp9 = c(7, 6, 3), sp10 = c(8, 7, 3), sp11 = c(2, 3, 2), sp12 = c(3, 3, 2))

  # Define state_var including the state variables for the states in the EDR and
  # the targets
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")],
                                target_var[, -c("traj", "state")]))

  # Define trajectories and states in state_var and the ID of the targets
  trajectories <- c(EDR_var$traj, target_var$traj)
  states <- as.integer(c(EDR_var$state, target_var$state))
  targets <- c("target1", "target2")

  # Define the function used to calculate state dissimilarities and its arguments
  d_function = "vegan::vegdist"
  d_args <- list(x = state_var, method = "bray")

  # Compute PETRA-EDR
  petra <- petra_edr(state_var = state_var,
                     trajectories = trajectories,
                     states = states,
                     targets = targets,
                     k = c(5L, 3L),
                     minPts = c(2L, 3L),
                     eps = c(NA, 0.6),
                     d_function = d_function,
                     d_args = d_args,
                     method = "mean",
                     w_function = c("exponential", NA),
                     alpha = c(3, NA))
}

```
