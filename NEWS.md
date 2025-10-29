# ecoregime 0.3.0

## Major changes

* `state_to_trajectory()`: `method = "projection"` considers the nearest states in `reference` to calculate dissimilarities between the `target_states` and `reference`. `method = "mixed"` returns the same values than `method = "projection"`.

## Minor changes and bug fixes

* `plot.RETRA()` correctly represents trajectories when states are not in order.

* `plot_edr()` does not print `[[i]]NULL` to the console.

* `retra_edr()`: fixed error when `!is.null(Dim)`.

* Fixed errors due to incompatibilities with ecotraj 1.1.1

# ecoregime 0.2.1

## Bugfixes

* Fixed errors due to incompatibilities with ecotraj 1.0.0

## New features

* `plot_edr()` represents EDR trajectories and states according to pre-defined 
colors or variables

# ecoregime 0.2.0

## New features

* `resistance()` calculates resistance to disturbances.

* `amplitude()` calculates the amplitude of disturbed trajectories.

* `recovery()` calculates the recovery towards reference trajectories.

* `net_change()` calculates the net change of disturbed trajectories.

* `state_to_trajectory()` defines the position of a state with respect to a trajectory.

* `EDR_data` now includes an abundance matrix for disturbed communities.

* New `vignette("Resilience")` which describes how ecological resilience can be assessed.

## Minor improvements

* `dEve()` error now does not refer to any reference trajectory.

* `define_retra()` returns error if certain characters are used in `trajectories`.

# ecoregime 0.1.3

* `dist_edr()` can be used when states are not in order.

# ecoregime 0.1.2

* Initial CRAN submission.
