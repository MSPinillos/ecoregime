dist_edr <- function(dtraj, edr, symmetrization, distance.type = "DDR"){

  distance.type <- match.arg(distance.type, c("DDR", "minDist", "maxDist"))

  # Convert dtraj into a matrix
  dtrajmat <- as.matrix(dtraj)

  # Check input data is correct
  if(length(edr) != nrow(dtrajmat)){
    stop(cat("'edr' needs to have a length equal to the number of rows and columns of 'dtraj'. Provide the EDR of each trajectory considered in 'dtraj'.", "\n"))
  }

  # Nb edr and trajectories/edr
  ID_edr <- unique(edr)
  Ntraj_edr <- table(edr)

  # Empty matrix to compile regime distances
  dreg <- matrix(0, ncol = length(ID_edr), nrow = length(ID_edr), dimnames = list(ID_edr, ID_edr))


  for(iR1 in seq_along(ID_edr)){                      # R1

    for(iR2 in seq_along(ID_edr)){                    # R2
      minD_T1R2 <- numeric()

      for(T1 in 1:Ntraj_edr[iR1]){        # for each T1i, find the minimum distance to R2; min{DSDSP(T1i, T21), ..., DSDSP(T1i, T2m)}
        iT1 <- sum(Ntraj_edr[1:iR1]) - Ntraj_edr[iR1] + T1                         # index of T1i in dtraj
        T2 <- (sum(Ntraj_edr[1:iR2]) - Ntraj_edr[iR2] + 1):sum(Ntraj_edr[1:iR2])  # indices of all T21,..., T2m in dtraj
        minD_T1R2 <- c(minD_T1R2, min(dtrajmat[iT1, T2]))            # min{DSDSP(T1i, T21), ..., DSDSP(T1i, T2m)}
      }
      if(distance.type == "DDR"){
        dreg[iR1, iR2] <- sum(minD_T1R2)/(Ntraj_edr[iR1])           # DDR(R1, R2)
      }
      if(distance.type == "minDist"){
        dreg[iR1, iR2] <- min(minD_T1R2)
      }
      if(distance.type == "maxDist"){
        dreg[iR1, iR2] <- max(minD_T1R2)
      }
    }
  }

  # Symmetrize DDR if required
  if(!is.null(symmetrization)){
    for(iR1 in seq_along(ID_edr)){
      for(iR2 in seq_along(ID_edr)){
        dreg[iR1, iR2] <- do.call(symmetrization, list(c(dreg[iR1, iR2], dreg[iR2, iR1])))
        dreg[iR2, iR1] <- dreg[iR1, iR2]
      }
    }
    dreg <- as.dist(dreg)
  }

  return(dreg)

}
