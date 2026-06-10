# Define representative trajectories from trajectory features

Generate an object of class `RETRA` from a data frame containing
trajectory states to define representative trajectories in Ecological
Dynamic Regimes (EDR).

## Usage

``` r
define_retra(data, d = NULL, trajectories = NULL, states = NULL, retra = NULL)
```

## Arguments

- data:

  A data frame of four columns indicating identifiers for the new
  representative trajectories, the individual trajectories or sites to
  which the states belong, the order of the states in the individual
  trajectories, and the identifier of the representative trajectory to
  which the states belong (only if `!is.null(retra)`). Alternatively,
  'data' can be a vector or a list of character vectors including the
  sequence of segments forming the new representative trajectory. See
  Details for further clarifications to define `data`.

- d:

  Either a symmetric matrix or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states of all trajectories in the
  EDR. If `NULL` (default), the length (`Length`) of the new
  representative trajectories and the distances between states of
  different trajectories or sites (`Link_distance`) are not calculated.

- trajectories:

  Only needed if `!is.null(d)`. Vector indicating the trajectory or site
  to which each state in `d` belongs.

- states:

  Only needed if `!is.null(d)`. Vector of integers indicating the order
  of the states in `d` for each trajectory.

- retra:

  Object of class `RETRA` returned from
  [`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md).
  If `NULL` (default), `minSegs` and `Seg_density` are not provided for
  the new representative trajectories.

## Value

An object of class `RETRA`, which is a list of length equal to the
number of representative trajectories defined. For each trajectory, the
following information is returned:

- `minSegs`:

  Value of the `minSegs` parameter used in
  [`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md).
  If `retra` is `NULL`, `minSegs` = `NA`.

- `Segments`:

  Vector of strings including the sequence of segments forming the
  representative trajectory. Each segment is identified by a string of
  the form `traj[st1-st2]`, where `traj` is the identifier of the
  original trajectory to which the segment belongs and `st1` and `st2`
  are identifiers of the initial and final states defining the segment.
  The same format `traj[st1-st2]` is maintained when only one state of
  an individual trajectory is considered (`st1` = `st2`). `traj`, `st1`,
  and `st2` are recycled from `data`.

- `Size`:

  Integer indicating the number of states forming the representative
  trajectory.

- `Length`:

  Numeric value indicating the length of the representative trajectory,
  calculated as the sum of the dissimilarities in `d` between every pair
  of consecutive states. If `d` is `NULL`, `Length` = `NA`.

- `Link_distance`:

  Data frame of two columns indicating artificial links between two
  segments (`Link`) and the dissimilarity between the connected states
  (`Distance`). When two representative segments are linked by a common
  state or by two consecutive states of the same trajectory, the link
  distance is zero or equal to the length of a real segment,
  respectively. In both cases, the link is not considered in the
  returned data frame. If `d` is `NULL`, `Link_distance` = `NA`.

- `Seg_density`:

  Data frame of two columns and one row for each representative segment.
  `Density` contains the number of segments in the EDR that is
  represented by each segment of the representative trajectory.
  `kdTree_depth` contains the depth of the k-d tree for each leaf
  represented by the corresponding segment. That is, the number of
  partitions of the ordination space until finding a region with
  `minSegs` segments or less. If `retra` is `NULL`, `Seg_density` =
  `NA`.

## Details

Each representative trajectory returned by the function
[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
corresponds to the longest sequence of representative segments that can
be linked according to the criteria defined in the RETRA-EDR algorithm
(Sánchez-Pinillos et al., 2023). One could be interested in splitting
the obtained trajectories, considering only a fraction of the returned
trajectories, or defining representative trajectories following
different criteria than those in RETRA-EDR. The function
`define_retra()` allows generating an object of class `RETRA` that can
be used in other functions of `ecoregime` (e.g.,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html)).

For that, it is necessary to provide information about the set of
segments or trajectory states that form the new representative
trajectory through the argument `data`:

- `data` can be defined as a **data frame** with as many rows as the
  number of states in all representative trajectories and the following
  columns:

  - `RT`:

    A string indicating the identifier of the new representative
    trajectories. Each identifier needs to appear as many times as the
    number of states forming each representative trajectory.

  - `RT_traj`:

    A vector indicating the individual trajectories in the EDR to which
    each state of the new representative trajectory belongs.

  - `RT_states`:

    A vector of integers indicating the identifier of the states forming
    the new representative trajectories. Each integer must refer to the
    order of the states in the individual trajectories of the EDR to
    which they belong.

  - `RT_retra`:

    Only if the new trajectories are defined from representative
    trajectories returned by
    [`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
    (when `!is.null(retra)`). A vector of strings indicating the
    representative trajectory in `retra` to which each state belongs.

- Alternatively, `data` can be defined as either a **vector** (if there
  is one representative trajectory) or a **list of character vectors**
  (with as many elements as the number of representative trajectories
  desired) containing the sequence of segments of the representative
  trajectories. In any case, each segment needs to be specified in the
  form `traj[st1-st2]`, where `traj` is the identifier of the original
  trajectory to which the segment belongs and `st1` and `st2` are
  identifiers of the initial and final states defining the segment. If
  only one state of an individual trajectory is considered to form the
  representative trajectory, the corresponding segment needs to be
  defined as `traj[st-st]`.

## See also

[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
for identifying representative trajectories in EDRs through RETRA-EDR.

[`summary()`](https://rdrr.io/r/base/summary.html) for summarizing the
characteristics of the representative trajectories.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for plotting
representative trajectories in an ordination space representing the
state space of the EDR.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
# Example 1 -----------------------------------------------------------------
# Define representative trajectories from the outputs of retra_edr().

# Identify representative trajectories using retra_edr()
d <- EDR_data$EDR1$state_dissim
trajectories <- EDR_data$EDR1$abundance$traj
states <- EDR_data$EDR1$abundance$state
old_retra <- retra_edr(d = d, trajectories = trajectories, states = states,
                       minSegs = 5)

# retra_edr() returns three representative trajectories
old_retra
#> $T1
#> $T1$minSegs
#> [1] 5
#> 
#> $T1$Segments
#> [1] "28[1-2]" "28[2-3]"
#> 
#> $T1$Size
#> [1] 3
#> 
#> $T1$Length
#> [1] 0.2669408
#> 
#> $T1$Link_distance
#> [1] NA
#> 
#> $T1$Seg_density
#>         Density kdTree_depth
#> 28[1-2]       7            4
#> 28[2-3]       6            5
#> 
#> 
#> $T2
#> $T2$minSegs
#> [1] 5
#> 
#> $T2$Segments
#> [1] "28[1-2]" "30[2-3]" "5[3-4]"  "15[1-2]" "4[2-3]"  "4[3-4]"  "1[1-2]" 
#> [8] "14[2-3]"
#> 
#> $T2$Size
#> [1] 15
#> 
#> $T2$Length
#> [1] 0.9270207
#> 
#> $T2$Link_distance
#>                Link    Distance
#> 1 28[1-2] - 30[2-3] 0.074626866
#> 2  30[2-3] - 5[3-4] 0.045685279
#> 3  5[3-4] - 15[1-2] 0.144278607
#> 4  15[1-2] - 4[2-3] 0.035175879
#> 5   4[3-4] - 1[1-2] 0.010000000
#> 6  1[1-2] - 14[2-3] 0.004926108
#> 
#> $T2$Seg_density
#>         Density kdTree_depth
#> 28[1-2]       7            4
#> 30[2-3]       6            5
#> 5[3-4]        6            5
#> 15[1-2]       6            6
#> 4[2-3]        7            7
#> 4[3-4]        9            5
#> 1[1-2]        9            3
#> 14[2-3]       7            6
#> 
#> 
#> $T3
#> $T3$minSegs
#> [1] 5
#> 
#> $T3$Segments
#> [1] "6[1-2]"  "5[3-4]"  "15[1-2]" "4[2-3]"  "4[3-4]"  "1[1-2]"  "14[2-3]"
#> 
#> $T3$Size
#> [1] 13
#> 
#> $T3$Length
#> [1] 0.6756596
#> 
#> $T3$Link_distance
#>               Link    Distance
#> 1  6[1-2] - 5[3-4] 0.020202020
#> 2 5[3-4] - 15[1-2] 0.144278607
#> 3 15[1-2] - 4[2-3] 0.035175879
#> 4  4[3-4] - 1[1-2] 0.010000000
#> 5 1[1-2] - 14[2-3] 0.004926108
#> 
#> $T3$Seg_density
#>         Density kdTree_depth
#> 6[1-2]        8            5
#> 5[3-4]        6            5
#> 15[1-2]       6            6
#> 4[2-3]        7            7
#> 4[3-4]        9            5
#> 1[1-2]        9            3
#> 14[2-3]       7            6
#> 
#> 
#> attr(,"class")
#> [1] "RETRA"

# Keep the last five segments of trajectories "T2" and "T3"
selected_segs <- old_retra$T2$Segments[4:length(old_retra$T2$Segments)]

# Identify the individual trajectories for each state...
selected_segs
#> [1] "15[1-2]" "4[2-3]"  "4[3-4]"  "1[1-2]"  "14[2-3]"
selected_traj <- rep(c(15, 4, 4, 1, 14), each = 2)

# ...and the states (in the same order as the representative trajectory).
selected_states <- c(1, 2, 2, 3, 3, 4, 1, 2, 2, 3)

# Generate the data frame with the format indicated in the documentation
df <- data.frame(RT = rep("A", length(selected_states)),
                 RT_traj = selected_traj,
                 RT_states = as.integer(selected_states),
                 RT_retra = rep("T2", length(selected_states)))

# Remove duplicates (trajectory 4, state 3)
df <- unique(df)

# Generate a RETRA object using define_retra()
new_retra <- define_retra(data = df,
                          d = d,
                          trajectories = trajectories,
                          states = states,
                          retra = old_retra)

# Example 2 -----------------------------------------------------------------
# Define representative trajectories from sequences of segments

# Select all segments in T1, split T2 into two new trajectories, and include
# a trajectory composed of states belonging to trajectories "5", "6", and "7"
data <- list(old_retra$T1$Segments,
             old_retra$T2$Segments[1:3],
             old_retra$T2$Segments[4:8],
             c("5[1-2]", "5[2-3]", "7[4-4]", "6[4-5]"))

# Generate a RETRA object using define_retra()
new_retra <- define_retra(data = data,
                          d = d,
                          trajectories = trajectories,
                          states = states,
                          retra = old_retra)

# Example 3 -----------------------------------------------------------------
# Define two representative trajectories from individual trajectories in EDR1.

# Define trajectory "A" from states in trajectories 3 and 4
data_A <- data.frame(RT = rep("A", 4),
                     RT_traj = c(3, 3, 4, 4),
                     RT_states = c(1:2, 4:5))

# Define trajectory "B" from states in trajectories 5, 6, and 7
data_B <- data.frame(RT = rep("B", 5),
                     RT_traj = c(5, 5, 7, 6, 6),
                     RT_states = c(1, 2, 4, 4, 5))

# Compile data for both trajectories in a data frame
df <- rbind(data_A, data_B)
df$RT_states <- as.integer(df$RT_states)

# Generate a RETRA object using define_retra()
new_retra <- define_retra(data = df, d = EDR_data$EDR1$state_dissim,
                          trajectories = EDR_data$EDR1$abundance$traj,
                          states = EDR_data$EDR1$abundance$state)

```
