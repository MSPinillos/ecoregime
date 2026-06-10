# Summarize representative trajectories

Summarize the properties of representative trajectories returned by
[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
or
[`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)

## Usage

``` r
# S3 method for class 'RETRA'
summary(object, ...)
```

## Arguments

- object:

  An object of class `RETRA`.

- ...:

  (not used)

## Value

Data frame with nine columns and one row for each representative
trajectory in `object`. The columns in the returned data frame contain
the following information:

- `ID`:

  Identifier of the representative trajectories.

- `Size`:

  Number of states forming each representative trajectory.

- `Length`:

  Sum of the dissimilarities in `d` between every pair of consecutive
  states forming the representative trajectories.

- `Avg_link`:

  Mean value of the dissimilarities between consecutive states of the
  representative trajectories that do not belong to the same ecological
  trajectory or site (i.e., artificial links).

- `Sum_link`:

  Sum of the dissimilarities between consecutive states of the
  representative trajectories that do not belong to the same ecological
  trajectory or site (i.e., artificial links).

- `Avg_density`:

  Mean value of the number of segments represented by each segment of
  the representative trajectory (excluding artificial links).

- `Max_density`:

  Maximum number of segments represented by at least one of the segments
  of the representative trajectory (excluding artificial links).

- `Avg_depth`:

  Mean value of the k-d tree depths, that is, the number of partitions
  of the ordination space until finding a region with `minSegs` segments
  or less.

- `Max_depth`:

  Maximum depth in the k-d tree, that is, the number of partitions of
  the ordination space until finding a region with `minSegs` segments or
  less.

## See also

[`retra_edr()`](https://mspinillos.github.io/ecoregime/reference/retra_edr.md)
for identifying representative trajectories in EDRs applying RETRA-EDR.

[`define_retra()`](https://mspinillos.github.io/ecoregime/reference/define_retra.md)
for generating an object of class `RETRA` from trajectory features.

## Examples

``` r
# Apply RETRA-EDR to identify representative trajectories
d = EDR_data$EDR1$state_dissim
trajectories = EDR_data$EDR1$abundance$traj
states = EDR_data$EDR1$abundance$state
RT <- retra_edr(d = d, trajectories = trajectories, states = states, minSegs = 5)

# Summarize the properties of the representative trajectories in a data frame
summary(RT)
#>    ID Size    Length   Avg_link  Sum_link Avg_density Max_density Avg_depth
#> T1 T1    3 0.2669408         NA        NA    6.500000           7  4.500000
#> T2 T2   15 0.9270207 0.05244879 0.3146927    7.125000           9  5.125000
#> T3 T3   13 0.6756596 0.04291652 0.2145826    7.428571           9  5.285714
#>    Max_depth
#> T1         5
#> T2         7
#> T3         7
```
