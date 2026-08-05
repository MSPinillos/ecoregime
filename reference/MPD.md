# Mean predicted deviation (MPD)

Metric to quantify the accuracy of predicted trajectories using
PETRA-EDR.

## Usage

``` r
MPD(x, pc_predicted = 1)
```

## Arguments

- x:

  Object of class `PETRA` returned by
  [`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md)
  using `return_args = T`.

- pc_predicted:

  Numeric value in the range 0-1, indicating the percentage of predicted
  states used to compute MPD.

## Value

`MPD()` returns a numeric value quantifying the average dissimilarity
between the target used to compute `x` and the trajectories predicted
from the precedent and following states included in `x`. Note that the
larger the value of MPD, the lower the accuracy of `x`.

## Details

MDP estimates the accuracy of predicted trajectories (`x`) returned by
[`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md).
MDP is based on the assumption that only one trajectory passes by each
target, and therefore, subsequent trajectories predicted from any point
of `x` should contain the original states in the target.

\\ MPD = \frac{\sum\_{i=1}^{n}\sum\_{j=1}^{m}d(x_i, \hat{Z}\_j)}{n m} \\

where \\d(x_i, \hat{Z}\_i)\\ is the dissimilarity between the initial or
final state \\x_i\\ of the target and the predicted trajectory
\\\hat{Z}\_j\\, defined as the minimum dissimilarity between \\x_i\\ and
the predicted states of \\\hat{Z}\_j\\; \\n\\ is the number of states of
the target; and \\m\\ is the number of states of the trajectory
predicted from the target (\\\hat{Z}\_i\\).

`pc_predicted` values lower than 1 may reduce the computation time.

## References

Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
Forecasting ecological trajectories from ecological dynamic regimes to
improve resilience analysis. *Methods in Ecology and Evolution*.
<doi:10.1111/2041-210x.70372>

## See also

[`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md)
for predicting ecological trajectories from EDRs.

[`plot.PETRA()`](https://mspinillos.github.io/ecoregime/reference/plot.PETRA.md)
for plotting predicted trajectories in an ordination space representing
the state space of the EDR.

## Author

Martina Sánchez-Pinillos

## Examples

``` r
if (requireNamespace("vegan", quietly = TRUE)) {
  # Define state_var including the state variables for the states in the EDR and
  # the target
  EDR_var <- EDR_data$EDR1$abundance
  target_var <- data.frame(sp1 = 30, sp2 = 10, sp3 = 13, sp4 = 12, sp5 = 5, sp6 = 2,
                           sp7 = 6, sp8 = 3, sp9 = 7, sp10 = 8, sp11 = 2, sp12 = 3)
  state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))

  # Compute PETRA-EDR
  petra <- petra_edr(state_var = state_var,
                     trajectories = c(EDR_var$traj, "target"),
                     states = as.integer(c(EDR_var$state, 0)),
                     targets = "target",
                     k = 5L,
                     minPts = 2L,
                     d_function = "vegan::vegdist",
                     d_args = list(x = state_var, method = "bray"),
                     return_args = TRUE,
                     direction = 2)

  # Calculate MPD
  MPD(x = petra)
}
#>     target 
#> 0.08324062 

```
