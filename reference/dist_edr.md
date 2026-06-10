# Dissimilarities between Ecological Dynamic Regimes

Generate a matrix containing dissimilarities between one or more pairs
of Ecological Dynamic Regimes (EDR). `dist_edr()` computes different
dissimilarity indices, all of them based on the dissimilarities between
the trajectories of two EDRs.

## Usage

``` r
dist_edr(
  d,
  d.type,
  trajectories = NULL,
  states = NULL,
  edr,
  metric = "dDR",
  symmetrize = NULL,
  ...
)
```

## Arguments

- d:

  Symmetric matrix or object of class
  [`dist`](https://rdrr.io/r/stats/dist.html) containing the
  dissimilarities between each pair of states of all trajectories in the
  EDR or the dissimilarities between each pair of trajectories.

- d.type:

  One of `"dStates"` (if `d` contains state dissimilarities) or
  `"dTraj"` (if `d` contains trajectory dissimilarities).

- trajectories:

  Only if `d.type` = `"dStates"`. Vector indicating the trajectory or
  site corresponding to each entry in `d`.

- states:

  Only if `d.type` = `"dStates"`. Vector of integers indicating the
  order of the states in `d` for each trajectory.

- edr:

  Vector indicating the EDR to which each trajectory/state in `d`
  belongs.

- metric:

  A string indicating the dissimilarity index to be used: `"dDR"`
  (default), `"minDist"`, `"maxDist"`.

- symmetrize:

  String naming the function to be called to symmetrize the resulting
  dissimilarity matrix (`"mean"`, `"min"`, `"max`, `"lower"`,
  `"upper"`). If `NULL` (default), the matrix is not symmetrized.

- ...:

  Only if `d.type` = `"dStates"`. Further arguments to calculate
  trajectory dissimilarities. See
  [`ecotraj::trajectoryDistances()`](https://emf-creaf.github.io/ecotraj/reference/trajectoryComparison.html).

## Value

Matrix including the dissimilarities between every pair of EDRs.

## Details

The implemented metrics are:

- `"dDR"`:

  \\ d\_{DR}(R_1, R_2) = \frac{1}{n} \sum\_{i=1}^{n} d\_{TR}(T\_{1i},
  R_2) \\

- `"minDist"`:

  \\ d\_{DRmin}(R_1, R_2) = \min\_{i=1}^{n} \\ d\_{TR}(T\_{1i}, R_2) \\
  \\

- `"maxDist"`:

  \\ d\_{DRmax}(R_1, R_2) = \max\_{i=1}^{n} \\ d\_{TR}(T\_{1i}, R_2) \\
  \\

where \\R_1\\ and \\R_2\\ are two EDRs composed of \\n\\ and \\m\\
ecological trajectories, respectively, and \\d\_{TR}(T\_{1i}, R_2)\\ is
the dissimilarity between the trajectory \\T\_{1i}\\ of \\R_1\\ and the
closest trajectory of \\R_2\\:

\\ d\_{TR}(T\_{1i}, R_2) = \min\\d_T(T\_{1i}, T\_{21}), ... ,
d_T(T\_{1i}, T\_{2m})\\ \\

The metrics calculated are not necessarily symmetric. That is,
\\d\_{DR}(R_1, R_2)\\ is not necessarily equal to \\d\_{DR}(R_2, R_1)\\.
It is possible to symmetrize the returned matrix by indicating the name
of the function to be used in `symmetrize`:

- `"mean"`:

  \\ d\_{DRsym} = \frac{d\_{DR}(R_1, R_2) + d\_{DR}(R_2, R_1)}{2} \\

- `"min"`:

  \\ d\_{DRsym} = \min\\d\_{DR}(R_1, R_2), d\_{DR}(R_2, R_1)\\ \\

- `"max"`:

  \\ d\_{DRsym} = \max\\d\_{DR}(R_1, R_2), d\_{DR}(R_2, R_1)\\ \\

- `"lower"`:

  The lower triangular part of the dissimilarity matrix is used.

- `"upper"`:

  The upper triangular part of the dissimilarity matrix is used.

## References

Sánchez-Pinillos, M., Kéfi, S., De Cáceres, M., Dakos, V. 2023.
Ecological Dynamic Regimes: Identification, characterization, and
comparison. *Ecological Monographs*. <doi:10.1002/ecm.1589>

## Author

Martina Sánchez-Pinillos

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE)) {
    # Load species abundances and compile in a data frame
    abun1 <- EDR_data$EDR1$abundance
    abun2 <- EDR_data$EDR2$abundance
    abun3 <- EDR_data$EDR3$abundance
    abun <- data.frame(rbind(abun1, abun2, abun3))

    # Define row names in abun to keep the reference of the EDR, trajectory, and
    # state
    row.names(abun) <- paste0(abun$EDR, "_", abun$traj, "_", abun$state)

    # Calculate dissimilarities between every pair of states
    # For example, Bray-Curtis index
    dStates <- vegan::vegdist(abun[, -c(1, 2, 3)], method = "bray")

    # Use the labels in dStates to define the trajectories to which each state
    # belongs
    id_traj <- vapply(strsplit(labels(dStates), "_"), function(x){
      paste0(x[1], "_", x[2])
    }, character(1))
    id_state <- vapply(strsplit(labels(dStates), "_"), function(x){
    as.integer(x[3])
    }, integer(1))
    id_edr <- vapply(strsplit(labels(dStates), "_"), function(x){
    paste0("EDR", x[1])
    }, character(1))

    # Calculate dissimilarities between every pair of trajectories
    dTraj <- ecotraj::trajectoryDistances(ecotraj::defineTrajectories(d = dStates, sites = id_traj,
                                                                      surveys = id_state),
                                          distance.type = "DSPD")

    # Use labels in dTraj to identify EDRs
    id_edr_traj <- vapply(strsplit(labels(dTraj), "_"), function(x){
    paste0("EDR", x[1])
    }, character(1))

    # Compute dissimilarities between EDRs:
    # 1) without symmetrizing the matrix and using state dissimilarities
    dEDR <- dist_edr(d = dStates, d.type = "dStates",
                     trajectories = id_traj, states = id_state, edr = id_edr,
                     metric = "dDR", symmetrize = NULL)

    # 2) symmetrizing by averaging elements on and below the diagonal and using
    # trajectory dissimilarities
    dEDR <- dist_edr(d = dTraj, d.type = "dTraj", edr = id_edr_traj,
                     metric = "dDR", symmetrize = "mean")
}
```
