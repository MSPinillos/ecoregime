# Representative trajectories in Ecological Dynamic Regimes (RETRA-EDR)

`retra_edr()` applies the algorithm RETRA-EDR (Sánchez-Pinillos et al.,
2023) to identify representative trajectories summarizing the main
dynamical patterns of an Ecological Dynamic Regime (EDR).

## Usage

``` r
retra_edr(
  d,
  trajectories,
  states,
  minSegs,
  dSegs = NULL,
  coordSegs = NULL,
  traj_Segs = NULL,
  state1_Segs = NULL,
  state2_Segs = NULL,
  Dim = NULL,
  eps = 0
)
```

## Arguments

- d:

  Either a symmetric matrix or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states of all trajectories in the
  EDR.

- trajectories:

  Vector indicating the trajectory or site to which each state in `d`
  belongs.

- states:

  Vector of integers indicating the order of the states in `d` for each
  trajectory.

- minSegs:

  Integer indicating the minimum number of segments in a region of the
  EDR represented by a segment of the representative trajectory.

- dSegs:

  Either a symmetric matrix or an object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between every pair of trajectory segments (see
  Details).

- coordSegs:

  Matrix containing the coordinates of trajectory segments (rows) in
  each axis (columns) of an ordination space (see Details).

- traj_Segs:

  Vector indicating the trajectory to which each segment in `dSeg`
  and/or `coordSegs` belongs. Only required if `dSegs` or `coordSegs`
  are not `NULL`.

- state1_Segs:

  Vector indicating the initial state of each segment in `dSegs` and/or
  `coordSegs` according to the values given in `states`. Only required
  if `dSegs` or `coordSegs` are not `NULL`.

- state2_Segs:

  Vector indicating the final state of each segment in `dSegs` and/or
  `coordSegs` according to the values given in `states`. Only required
  if `dSegs` or `coordSegs` are not `NULL`.

- Dim:

  Optional integer indicating the number of axes considered to partition
  the segment space and generate a k-d tree. By default (`Dim = NULL`),
  all axes are considered.

- eps:

  Numeric value indicating the minimum length in the axes of the segment
  space to be partitioned when the k-d tree is generated. If `eps = 0`
  (default), partitions are made regardless of the size.

## Value

The function `retra_edr()` returns an object of class `RETRA`, which is
a list of length equal to the number of representative trajectories
identified. For each trajectory, the following information is returned:

- `minSegs`:

  Value of the `minSegs` parameter.

- `Segments`:

  Vector of strings including the sequence of segments forming the
  representative trajectory. Each segment is identified by a string of
  the form `traj[st1-st2]`, where `traj` is the identifier of the
  original trajectory to which the segment belongs and `st1` and `st2`
  are identifiers of the initial and final states defining the segment.

- `Size`:

  Numeric value indicating the number of states forming the
  representative trajectory.

- `Length`:

  Numeric value indicating the length of the representative trajectory,
  calculated as the sum of the dissimilarities in `d` between every pair
  of consecutive states.

- `Link_distance`:

  Data frame of two columns indicating artificial links between
  representative segments (`Link`) and the dissimilarity between the
  connected states (`Distance`). When two representative segments are
  linked by a common state or by two consecutive states of the same
  trajectory, the link distance is zero or equal to the length of a real
  segment, respectively. In both cases, the link is not considered in
  the returned data frame.

- `Seg_density`:

  Data frame of two columns and one row for each representative segment.
  `Density` contains the number of segments in the EDR that is
  represented by each segment of the representative trajectory.
  `kdTree_depth` contains the depth of the k-d tree for each leaf
  represented by the corresponding segment. That is, the number of
  partitions of the ordination space until finding a region with
  `minSegs` segments or less.

## Details

The algorithm RETRA-EDR is based on a partition-and-group approach by
which it identifies regions densely crossed by ecological trajectories
in an EDR, selects a representative segment in each dense region, and
joins the representative segments by a set of artificial `Links` to
generate a network of representative trajectories. For that, RETRA-EDR
splits the trajectories of the EDR into segments and uses an ordination
space generated from a matrix containing the dissimilarities between
trajectory segments. Dense regions are identified by applying a k-d tree
to the ordination space.

By default, RETRA-EDR calculates segment dissimilarities following the
approach by De Cáceres et al. (2019) and applies metric multidimensional
scaling (mMDS, Borg and Groenen, 2005) to generate the ordination space.
It is possible to use other dissimilarity metrics and/or ordination
methods and reduce the computational time by indicating the
dissimilarity matrix and the coordinates of the segments in the
ordination space through the arguments `dSegs` and `coordSegs`,
respectively.

- If `!is.null(dSegs)` and `is.null(coordSegs)`, RETRA-EDR is computed
  by applying mMDS to `dSegs`.

- If `!is.null(dSegs)` and `!is.null(coordSegs)`, RETRA-EDR is directly
  computed from the coordinates provided in `coordSegs` and
  representative segments are identified using `dSegs`. `coordSegs`
  should be calculated by the user from `dSegs`.

- If `is.null(dSegs)` and `!is.null(coordSegs)` (not recommended),
  RETRA-EDR is directly computed from the coordinates provided in
  `coordSegs`. As `dSegs` is not provided, `retra_edr()` assumes that
  the ordination space is metric and identifies representative segments
  using the Euclidean distance.

## References

Borg, I., & Groenen, P. J. F. (2005). Modern Multidimensional Scaling
(2nd ed.). Springer.

De Cáceres, M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit
R & Hubbell S. (2019). Trajectory analysis in community ecology.
Ecological Monographs.

Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023.
Ecological Dynamic Regimes: Identification, characterization, and
comparison. *Ecological Monographs*. <doi:10.1002/ecm.1589>

## See also

[`summary()`](https://rdrr.io/r/base/summary.html) for summarizing the
characteristics of the representative trajectories.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) for plotting
representative trajectories in an ordination space representing the
state space of the EDR.

[`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)
for defining representative trajectories from a subset of segments or
trajectory features.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
# Example 1 -----------------------------------------------------------------
# Identify representative trajectories from state dissimilarities

# State dissimilarities (Bray-Curtis) from species abundances
abundance <- data.frame(EDR_data$EDR1$abundance)
d <- EDR_data$EDR1$state_dissim
# (d is equivalent to vegan::vegdist(abundance[, -c(1:3)]))

# Identify the trajectory (or site) and states in d
trajectories <- abundance$traj
states <- as.integer(abundance$state)

# Compute RETRA-EDR
RT1 <- retra_edr(d = d, trajectories = trajectories, states = states,
                 minSegs = 5)

# Example 2 -----------------------------------------------------------------
# Identify representative trajectories from segment dissimilarities

# Calculate segment dissimilarities using the Hausdorff distance
dSegs <- ecotraj::segmentDistances(ecotraj::defineTrajectories(d = d, sites = trajectories,
                                                               surveys = states),
                                   distance.type = "Hausdorff")
dSegs <- dSegs$Dseg

# Identify the trajectory (or site) and states in dSegs:
# Split the labels of dSegs (traj[st1-st2]) into traj, st1, and st2
seg_components <- strsplit(gsub("\\]", "", gsub("\\[", "-", labels(dSegs))), "-")
traj_Segs <- sapply(seg_components, "[", 1)
state1_Segs <- as.integer(sapply(seg_components, "[", 2))
state2_Segs <- as.integer(sapply(seg_components, "[", 3))

# Compute RETRA-EDR
RT2 <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5,
                dSegs = dSegs, traj_Segs = traj_Segs,
                state1_Segs = state1_Segs, state2_Segs = state2_Segs)

```
