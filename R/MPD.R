#' Mean predicted deviation (MPD)
#'
#' @description
#' Metric to quantify the accuracy of predicted trajectories using PETRA-EDR.
#'
#' @param x Object of class `PETRA` returned by `petra_edr()` using `return_args = T`.
#' @param pc_predicted Numeric value in the range 0-1, indicating the percentage
#' of predicted states used to compute MPD.
#'
#' @details
#' MDP estimates the accuracy of predicted trajectories (`x`) returned by `petra_edr()`.
#' MDP is based on the assumption that only one trajectory passes by each target,
#' and therefore, subsequent trajectories predicted from any point of `x` should
#' contain the original states in the target.
#'
#' \eqn{
#' MPD = \frac{\sum_{i=1}^{n}\sum_{j=1}^{m}d(x_i, \hat{Z}_j)}{n m}
#' }
#'
#' where \eqn{d(x_i, \hat{Z}_i)} is the dissimilarity between the initial or final
#' state \eqn{x_i} of the target and the predicted trajectory \eqn{\hat{Z}_j},
#' defined as the minimum dissimilarity between \eqn{x_i} and the predicted states
#' of \eqn{\hat{Z}_j}; \eqn{n} is the number of states of the target; and \eqn{m}
#' is the number of states of the trajectory predicted from the target (\eqn{\hat{Z}_i}).
#'
#' `pc_predicted` values lower than 1 may reduce the computation time.
#'
#' @returns
#' `MPD()` returns a numeric value quantifying the average dissimilarity between
#' the target used to compute `x` and the trajectories predicted from the precedent
#' and following states included in `x`. Note that the larger the value of MPD, the
#' lower the accuracy of `x`.
#'
#' @author Martina Sánchez-Pinillos
#'
#' @references
#' Sánchez-Pinillos M., Fortin, M-J., Messier, C., Kneeshaw, D. 2026.
#' Forecasting ecological trajectories from ecological dynamic regimes to improve
#' resilience analysis. *Methods in Ecology and Evolution*. <doi:10.1111/2041-210x.70372>
#'
#' @seealso
#' [`petra_edr()`] for predicting ecological trajectories from EDRs.
#'
#' [`plot.PETRA()`] for plotting predicted trajectories in an ordination space
#' representing the state space of the EDR.
#'
#' @export
#'
#' @examples
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   # Define state_var including the state variables for the states in the EDR and
#'   # the target
#'   EDR_var <- EDR_data$EDR1$abundance
#'   target_var <- data.frame(sp1 = 30, sp2 = 10, sp3 = 13, sp4 = 12, sp5 = 5, sp6 = 2,
#'                            sp7 = 6, sp8 = 3, sp9 = 7, sp10 = 8, sp11 = 2, sp12 = 3)
#'   state_var <- data.frame(rbind(EDR_var[, -c("EDR", "traj", "state")], target_var))
#'
#'   # Compute PETRA-EDR
#'   petra <- petra_edr(state_var = state_var,
#'                      trajectories = c(EDR_var$traj, "target"),
#'                      states = as.integer(c(EDR_var$state, 0)),
#'                      targets = "target",
#'                      k = 5L,
#'                      minPts = 2L,
#'                      d_function = "vegan::vegdist",
#'                      d_args = list(x = state_var, method = "bray"),
#'                      return_args = TRUE,
#'                      direction = 2)
#'
#'   # Calculate MPD
#'   MPD(x = petra)
#' }
#'
#'
MPD <- function(x, pc_predicted = 1) {

  # Check the format of x
  if (!inherits(x, "PETRA")) {
    stop("'x' must be an object of class 'PETRA'.")
  }
  if (!"arguments" %in% names(x)) {
    stop("'x' must contain 'arguments'. Set `return_arguments = T` when computing 'petra_edr()' to get 'x'.")
  }

  ## DEFINE VARIABLES ----------------------------------------------------------

  # Variables used to compute x
  state_var <- x$arguments$d_args[[x$arguments$args_state_var]]
  trajectories <- x$arguments$trajectories
  states <- x$arguments$states
  targets <- x$arguments$targets

  ##----------------------------------------------------------------------------

  mpd <- vapply(setNames(targets, targets), function(itarget){
    # ID of the target in the original data
    ID_target <- which(trajectories %in% itarget)

    # ID of the target in x
    ID_target_x <- which(x$trajectories %in% itarget)

    # ID of a random sample of the predicted states in x
    x_trajst <- paste0(x$trajectories[ID_target_x], "_", x$states[ID_target_x])
    trajst <- paste0(trajectories[ID_target], "_", states[ID_target])
    x_trajst <- c(sample(x_trajst[!x_trajst %in% trajst],
                         round(pc_predicted*(length(x_trajst) - length(trajst)))),
                  trajst)
    x_trajst <- x_trajst[order(x_trajst)]
    ID_predicted <- ID_target_x[!(x_trajst %in% trajst)]

    if (length(ID_predicted) > 0) {
      # Remove the targets from the original data and add predicted states
      if (inherits(state_var, "list")) {
        new_state_var <- c(state_var[-ID_target], x$state_var[ID_predicted])
        class(new_state_var) <- class(state_var)
      } else {
        new_state_var <- rbind(state_var[-ID_target, ], x$state_var[ID_predicted, ])
      }
      new_trajectories <- c(trajectories[-ID_target],
                            paste0(x$trajectories[ID_predicted], "_", x$states[ID_predicted]))
      new_states <- as.integer(c(states[-ID_target], x$states[ID_predicted]))
      new_targets <- paste0(x$trajectories[ID_predicted], "_", x$states[ID_predicted])
      if (!is.null(x$arguments$w)) {
        if (inherits(x$arguments$w, "list")) {
          new_w <- c(x$arguments$w[[which(targets == itarget)]][-ID_target], rep(0, length(ID_predicted)))
        } else {
          new_w <- c(x$arguments$w[-ID_target], rep(0, length(ID_predicted)))
        }
      } else {
        new_w <- NULL
      }

      # Identify the argument in d_function corresponding to state_var
      istate_var <- x$arguments$args_state_var

      # Dissimilarity between the predicted states and each state of the EDR
      d_predicted <- lapply(new_targets, function(jtarget){
        d_function <- x$arguments$d_function
        d_args <- x$arguments$d_args

        if (inherits(new_state_var, "list")) {
          vapply(1:length(new_state_var), function(istate){
            d_args[[istate_var]] <- c(new_state_var[istate], new_state_var[new_trajectories == jtarget])
            class(d_args[[istate_var]]) <- class(state_var)
            do.call(eval(parse(text = d_function)), args = d_args)
          }, numeric(1))
        } else {
          vapply(1:nrow(new_state_var), function(istate){
            d_args[[istate_var]] <- rbind(new_state_var[istate, ], new_state_var[new_trajectories == jtarget, ])
            do.call(eval(parse(text = d_function)), args = d_args)
          }, numeric(1))
        }
      })

      # Dissimilarity matrix
      new_d <- as.matrix(x$arguments$d)[-ID_target, -ID_target]
      for (id in d_predicted) {
        new_d <- rbind(new_d, id[1:ncol(new_d)])
      }
      for (id in d_predicted) {
        new_d <- cbind(new_d, id)
      }

      # Identify precedent/following predicted states in new_targets
      ID_pre_target <- which(x$states[ID_predicted] < min(states[ID_target]))
      ID_post_target <- which(x$states[ID_predicted] > max(states[ID_target]))

      # Define d_args
      new_d_args <- x$arguments$d_args
      new_d_args[[istate_var]] <- new_state_var

      # Apply PETRA-EDR to the preceding predicted states
      if (length(ID_pre_target) > 0) {
        y_pre = petra_edr(state_var = new_state_var,
                          trajectories = new_trajectories,
                          states = new_states,
                          targets = new_targets[ID_pre_target],
                          k = x$arguments$k[which(targets == itarget)],
                          eps = x$arguments$eps[which(targets == itarget)],
                          minPts = x$arguments$minPts[which(targets == itarget)],
                          d_function = x$arguments$d_function,
                          d_args = new_d_args,
                          d = new_d,
                          method = x$arguments$method,
                          w = new_w,
                          direction = 1)
      } else {
        y_pre <- NULL
      }

      # Apply PETRA-EDR to the following predicted states
      if (length(ID_post_target) > 0) {
        y_post = petra_edr(state_var = new_state_var,
                           trajectories = new_trajectories,
                           states = new_states,
                           targets = new_targets[ID_post_target],
                           k = x$arguments$k[which(targets == itarget)],
                           eps = x$arguments$eps[which(targets == itarget)],
                           minPts = x$arguments$minPts[which(targets == itarget)],
                           d_function = x$arguments$d_function,
                           d_args = new_d_args,
                           d = new_d,
                           method = x$arguments$method,
                           w = new_w,
                           direction = -1)
      } else {
        y_post <- NULL
      }


      # Calculate the minimum distance between the target states and the trajectories
      # predicted from the predicted states in x
      min_dist <- sapply(ID_target, function(itarget_st){
        # Compile data of the preceding and following states
        if (inherits(state_var, "list")) {
          state_var_ist <- c(state_var[itarget_st],
                             y_pre$state_var,
                             y_post$state_var)
          class(state_var_ist) <- class(state_var)
        } else {
          state_var_ist <- rbind(state_var[itarget_st, ],
                                 y_pre$state_var,
                                 y_post$state_var)
        }
        trajectories_ist <- c(trajectories[itarget_st],
                              y_pre$trajectories,
                              y_post$trajectories)
        states_ist <- c(states[itarget_st],
                        y_pre$states,
                        y_post$states)

        # Define d_args
        args_target_pred <- x$arguments$d_args

        sapply(unique(trajectories_ist)[-1], function(iTpre){
          # Define state_var using the target state and each of the predicted trajectories
          traj_target_pred <- trajectories_ist[trajectories_ist %in% c(trajectories_ist[1], iTpre)]
          st_target_pred <- states_ist[trajectories_ist %in% c(trajectories_ist[1], iTpre)]
          if (inherits(state_var_ist, "list")) {
            args_target_pred[[istate_var]] <- state_var_ist[which(trajectories_ist %in% traj_target_pred)]
            class(args_target_pred[[istate_var]]) <- class(state_var_ist)
          } else {
            args_target_pred[[istate_var]] <- state_var_ist[which(trajectories_ist %in% traj_target_pred), ]
          }

          # Compute dissimilarities
          d_target_pred <- as.matrix(do.call(eval(parse(text = x$arguments$d_function)),
                                             args = args_target_pred))

          if (length(traj_target_pred[-1]) > 1) {
            state_to_trajectory(d = d_target_pred,
                                trajectories = traj_target_pred,
                                states = as.integer(st_target_pred),
                                target_states = 1L,
                                reference = iTpre,
                                method = "nearest_state")$distance
          } else {
            d_target_pred[1, 2]
          }
        })
      })

      # MPD: average dissimilarities between the target states and the predicted
      # trajectories
      mpd_target <- mean(min_dist)

    } else {
      # If it was not possible to predict trajectories from the predicted states,
      # return NA
      mpd_target <- NA
    }
    return(mpd_target)
  }, numeric(1))

  return(mpd)

}

