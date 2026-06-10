# Metrics of Ecological Dynamic Regimes

Set of metrics to analyze the distribution and variability of
trajectories in Ecological Dynamic Regimes (EDR), including dynamic
dispersion (dDis), dynamic beta diversity (dBD), and dynamic evenness
(dEve).

## Usage

``` r
dDis(
  d,
  d.type,
  trajectories,
  states = NULL,
  reference,
  w.type = "none",
  w.values,
  ...
)

dBD(d, d.type, trajectories, states = NULL, ...)

dEve(d, d.type, trajectories, states = NULL, w.type = "none", w.values, ...)
```

## Arguments

- d:

  Symmetric matrix or object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states of all trajectories in the
  EDR or the dissimilarities between each pair of trajectories. To
  compute dDis, `d` needs to include the dissimilarities between all
  states/trajectories and the states/trajectory of reference.

- d.type:

  One of `"dStates"` (if `d` contains state dissimilarities) or
  `"dTraj"` (if `d` contains trajectory dissimilarities).

- trajectories:

  Vector indicating the trajectory or site corresponding to each entry
  in `d`.

- states:

  Only if `d.type` = `"dStates"`. Vector of integers indicating the
  order of the states in `d` for each trajectory.

- reference:

  Vector of the same class as `trajectories` and length equal to one,
  indicating the reference trajectory to compute dDis.

- w.type:

  Method used to weight individual trajectories:

  - `"none"`: All trajectories are considered equally relevant
    (default).

  - `"length"`: Trajectories are weighted by their length, calculated as
    the sum of the dissimilarities between every pair of consecutive
    states. `d` must contain dissimilarities between trajectory states
    and `d.type` = `"dStates"`.

  - `"size"`: Trajectories are weighted by their size, calculated as the
    number of states forming the trajectory. `d` must contain
    dissimilarities between trajectory states and `d.type` =
    `"dStates"`.

  - `"precomputed"`: Trajectories weighted according to different
    criteria.

- w.values:

  Only if `w.type` = `"precomputed"`. Numeric vector of length equal to
  the number of trajectories containing the weight of each trajectory.

- ...:

  Only if `d.type` = `"dStates"`. Further arguments to calculate
  trajectory dissimilarities. See
  [`ecotraj::trajectoryDistances()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.html).

## Value

- `dDis()` returns the value of dynamic dispersion for a given
  trajectory taken as a reference.

- `dBD()` returns the value of dynamic beta diversity.

- `dEve()` returns the value of dynamic evenness.

## Details

**Dynamic dispersion (`dDis()`)**

*dDis* is calculated as the average dissimilarity between each
trajectory in an EDR and a target trajectory taken as reference
(Sánchez-Pinillos et al., 2023).

\\ dDis = \frac{\sum\_{i=1}^{m}d\_{i\alpha}}{m} \\

where \\d\_{i\alpha}\\ is the dissimilarity between trajectory \\i\\ and
the trajectory of reference \\\alpha\\, and \\m\\ is the number of
trajectories.

Alternatively, it is possible to calculate a weighted mean of the
dissimilarities by assigning a weight to each trajectory.

\\ dDis =
\frac{\sum\_{i=1}^{m}w\_{i}d\_{i\alpha}}{\sum\_{i=1}^{m}w\_{i}} \\

where \\w\_{i}\\ is the weight assigned to trajectory \\i\\.

**Dynamic beta diversity (`dBD()`)**

*dBD* quantifies the overall variation of the trajectories in an EDR and
is equivalent to the average distance to the centroid of the EDR (De
Cáceres et al., 2019).

\\ dBD = \frac{\sum\_{i=1}^{m-1}\sum\_{j=i+1}^{m}d\_{ij}^{2}}{m(m-1)} \\

**Dynamic evenness (`dEve()`)**

*dEve* quantifies the regularity with which an EDR is filled by the
individual trajectories (Sánchez-Pinillos et al., 2023).

\\ dEve =
\frac{\sum\_{l=1}^{m-1}\min(\frac{d\_{ij}}{\sum\_{l=1}^{m-1}d\_{ij}},
\frac{1}{m-1}) - \frac{1}{m-1}}{1-\frac{1}{1-1}} \\

where \\d\_{ij}\\ is the dissimilarity between trajectories \\i\\ and
\\j\\ linked in a minimum spanning tree by the link \\l\\.

Optionally, it is possible to weight the trajectories of the EDR. In
that case, *dEve* becomes analogous to the functional evenness index
proposed by Villéger et al. (2008).

\\ dEve\_{w} =
\frac{\sum\_{l=1}^{m-1}\min(\frac{EW\_{ij}}{\sum\_{l=1}^{m-1}EW\_{ij}},
\frac{1}{m-1}) - \frac{1}{m-1}}{1-\frac{1}{1-1}} \\

where \\EW\_{ij}\\ is the weighted evenness:

\\ EW\_{ij} = \frac{d\_{ij}}{w_i + w_j} \\

## References

De Cáceres, M, Coll L, Legendre P, Allen RB, Wiser SK, Fortin MJ, Condit
R & Hubbell S. (2019). Trajectory analysis in community ecology.
Ecological Monographs.

Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023.
Ecological Dynamic Regimes: Identification, characterization, and
comparison. *Ecological Monographs*. <doi:10.1002/ecm.1589>

Villéger, S., Mason, N.W.H., Mouillot, D. (2008) New multidimensional
functional diversity indices for a multifaced framework in functional
ecology. Ecology.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
# Data to compute dDis, dBD, and dEve
dStates <- EDR_data$EDR1$state_dissim
dTraj <- EDR_data$EDR1$traj_dissim
trajectories <- paste0("T", EDR_data$EDR1$abundance$traj)
states <- EDR_data$EDR1$abundance$state

# Dynamic dispersion taking the first trajectory as reference
dDis(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories),
         reference = "T1")
#> dDis (ref. T1) 
#>       0.267622 

# Dynamic dispersion weighting trajectories by their length
dDis(d = dStates, d.type = "dStates", trajectories = trajectories, states = states,
         reference = "T1", w.type = "length")
#> dDis (ref. T1) 
#>      0.3548963 

# Dynamic beta diversity using trajectory dissimilarities
dBD(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories))
#>        dBD 
#> 0.03969095 

# Dynamic evenness
dEve(d = dStates, d.type = "dStates", trajectories = trajectories, states = states)
#>      dEve 
#> 0.7100773 

# Dynamic evenness considering that the 10 first trajectories are three times
# more relevant than the rest
w.values <- c(rep(3, 10), rep(1, length(unique(trajectories))-10))
dEve(d = dTraj, d.type = "dTraj", trajectories = unique(trajectories),
         w.type = "precomputed", w.values = w.values)
#>      dEve 
#> 0.7392367 
```
