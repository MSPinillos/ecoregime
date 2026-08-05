# Changelog

## ecoregime 0.4.1

CRAN release: 2026-08-04

### Minor changes and bug fixes

- Reduced computing time for examples.

## ecoregime 0.4.0

### New features

- [`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md)predicts
  longer trajectories of a target in ecological dynamic regimes.

- [`MPD()`](https://mspinillos.github.io/ecoregime/reference/MPD.md)
  quantifies the accuracy of predicted trajectories returned by
  [`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md).

- [`plot.PETRA()`](https://mspinillos.github.io/ecoregime/reference/plot.PETRA.md)
  represents predicted trajectories returned by
  [`petra_edr()`](https://mspinillos.github.io/ecoregime/reference/petra_edr.md)
  and the trajectories of the reference EDR in a state space.

- `vignette("Predicted_trajectories")`: describes how to predict
  ecological trajectories from EDRs.

### Improvements

- [`amplitude()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the amplitude of disturbed trajectories from `RETRA` and
  `PETRA` objects.

- [`recovery()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the recovery towards reference trajectories from `RETRA`
  and `PETRA` objects.

- [`net_change()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the net change of disturbed trajectories from `RETRA` and
  `PETRA` objects.

## ecoregime 0.3.1

CRAN release: 2026-06-07

### Minor changes and bug fixes

- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`plot_edr()`](https://mspinillos.github.io/ecoregime/reference/plot_edr.md):
  axes are not constrained to “Axis X”.

- [`vignette("EDR_framework")`](https://mspinillos.github.io/ecoregime/articles/EDR_framework.md):
  new plots showing ecological succession.

- [`vignette("Resilience")`](https://mspinillos.github.io/ecoregime/articles/Resilience.md):
  new Section 2.

- Packages in Suggests are used conditionally in examples, tests, and
  vignettes.

## ecoregime 0.3.0

CRAN release: 2025-10-29

### Major changes

- [`state_to_trajectory()`](https://mspinillos.github.io/ecoregime/reference/state_to_trajectory.md):
  `method = "projection"` considers the nearest states in `reference` to
  calculate dissimilarities between the `target_states` and `reference`.
  `method = "mixed"` returns the same values than
  `method = "projection"`.

### Minor changes and bug fixes

- [`plot.RETRA()`](https://mspinillos.github.io/ecoregime/reference/plot.RETRA.md)
  correctly represents trajectories when states are not in order.

- [`plot_edr()`](https://mspinillos.github.io/ecoregime/reference/plot_edr.md)
  does not print `[[i]]NULL` to the console.

- [`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md):
  fixed error when `!is.null(Dim)`.

- Fixed errors due to incompatibilities with ecotraj 1.1.1

## ecoregime 0.2.1

CRAN release: 2025-05-01

### Bugfixes

- Fixed errors due to incompatibilities with ecotraj 1.0.0

### New features

- [`plot_edr()`](https://mspinillos.github.io/ecoregime/reference/plot_edr.md)
  represents EDR trajectories and states according to pre-defined colors
  or variables

## ecoregime 0.2.0

CRAN release: 2024-04-17

### New features

- [`resistance()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates resistance to disturbances.

- [`amplitude()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the amplitude of disturbed trajectories.

- [`recovery()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the recovery towards reference trajectories.

- [`net_change()`](https://mspinillos.github.io/ecoregime/reference/resilience_metrics.md)
  calculates the net change of disturbed trajectories.

- [`state_to_trajectory()`](https://mspinillos.github.io/ecoregime/reference/state_to_trajectory.md)
  defines the position of a state with respect to a trajectory.

- `EDR_data` now includes an abundance matrix for disturbed communities.

- New
  [`vignette("Resilience")`](https://mspinillos.github.io/ecoregime/articles/Resilience.md)
  which describes how ecological resilience can be assessed.

### Minor improvements

- [`dEve()`](https://mspinillos.github.io/ecoregime/reference/EDR_metrics.md)
  error now does not refer to any reference trajectory.

- [`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)
  returns error if certain characters are used in `trajectories`.

## ecoregime 0.1.3

CRAN release: 2023-09-10

- [`dist_edr()`](https://mspinillos.github.io/ecoregime/reference/dist_edr.md)
  can be used when states are not in order.

## ecoregime 0.1.2

CRAN release: 2023-08-11

- Initial CRAN submission.
