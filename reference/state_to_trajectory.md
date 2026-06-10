# Position of a state with respect to a trajectory

Define the position of a state with respect to a reference trajectory
based on its distance from the trajectory and the length and direction
of the trajectory.

## Usage

``` r
state_to_trajectory(
  d,
  trajectories,
  states,
  target_states,
  reference,
  method,
  coordStates = NULL
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
  trajectory (assign 1 if the state does not belong to any trajectory).

- target_states:

  Vector of integers indicating the indices in `trajectories` and
  `states` of the ecological states for which their relative position
  will be calculated.

- reference:

  Vector of the same class of `trajectories` or object of class `RETRA`
  indicating the reference trajectory to calculate the relative position
  of the `target_states`.

- method:

  Method to calculate the distance and relative position of the
  `target_states` and the `reference`. One of `"nearest_state"` or
  `"projection"` (see Details).

- coordStates:

  Matrix containing the coordinates of each state (rows) and axis
  (columns) of a metric ordination space (see Details).

## Value

The function `state_to_trajectory()` returns a data frame of four
columns including the `distance` and `relative_position` between the
`target_state` and the `reference`.

- Depending on the `method`, `distance` is calculated as the
  dissimilarity between the `target_states` and their respective nearest
  state in the `reference` or the dissimilarity to their projections
  onto the `reference`.

- The `relative_position` is a value that ranges between 0 (if the
  nearest state or projected point coincides with the first `reference`
  state) and 1 (if the nearest state or projected point coincides with
  the last `reference` state).

## Details

`state_to_trajectory()` can calculate the distance and relative position
of one or more `target_states` relative to a `reference` trajectory by
two different methods:

- `"nearest_state"` returns the dissimilarity of the `target_states` to
  the nearest state of the `reference` trajectory (`distance`) and
  calculates the relative position of the nearest state within the
  `reference`.

- `"projection"` returns the dissimilarity of the `target_states` to
  their projection onto the `reference` trajectory and calculates the
  relative position of the projected state within the `reference`. When
  the `target_states` cannot be projected onto any of the segments
  forming the `reference` and in cases in which the dissimilarity to
  nearest state of the `reference` is smaller than the dissimilarity to
  the projected state, `state_to_trajectory()` uses the nearest state in
  the `reference` to compute `distance` and `relative_position`. This
  method requires `d` to be metric (i.e. to satisfy the triangle
  inequality). If `d` is not metric, `state_to_trajectory()` calculates
  the Euclidean distance within a transformed space generated through
  multidimensional scaling (Borg and Groenen, 2005). To use the state
  coordinates in a different metric space, use the `coordStates`
  argument.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
# State dissimilarities, trajectories, and states
d <- EDR_data$EDR3$state_dissim
# (d is Equivalent to vegan::vegdist(EDR_data$EDR3$abundance[, -c(1:3)]))
trajectories <- EDR_data$EDR3$abundance$traj
states <- EDR_data$EDR3$abundance$state

# Calculate the representative trajectories of an EDR to be used as reference
RT <- retra_edr(d = d,
               trajectories = trajectories,
               states = states,
               minSegs = 10)

# Define the target states
target_states <- as.integer(c(1, 16, 55))

# Calculate the position of the target states with respect to  the representative
# trajectories of an EDR
state_to_trajectory(d = d, trajectories = trajectories,
                    states = states,
                    target_states = target_states,
                    reference = RT,
                    method = "nearest_state")
#>   target_state reference   distance relative_position
#> 1            1        T1 0.04522613         0.0000000
#> 2            1        T2 0.75879397         0.3418326
#> 3           16        T1 0.07537688         0.7195150
#> 4           16        T2 0.07537688         0.3418326
#> 5           55        T1 0.04477612         0.2328098
#> 6           55        T2 0.58000000         0.3418326
```
