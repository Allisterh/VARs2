# ************************************************************************
# Bootstrap IVAR ----
# ************************************************************************
bootstrapIVAR <- function(IVAR, state, nboot, alpha = 90, method = "residual", fast = TRUE) {
  # unpacking
  data <- IVAR$data
  c_case <- IVAR$c_case
  t <- IVAR$t
  ub <- IVAR$u
  n <- IVAR$n
  u <- IVAR$u
  A <- IVAR$A
  p <- IVAR$p
  shockpos <- state$shockpos
  statepos <- state$statepos
  h <- IVAR$GIRF1$h
  impulse <- IVAR$GIRF1$impulse

  # define upper and lower thresholds
  lo <- (100 - alpha) / 2
  up <- 100 - lo

  GIRF1b <- array(NA, dim = c(h, n, nboot))
  GIRF2b <- array(NA, dim = c(h, n, nboot))

  nexplos <- 0
  nattempts <- 1
  print("Progress")
  bb <- 1
  pb <- txtProgressBar(min = 0, max = nboot, initial = 0, style = 3)
  while (bb <= nboot) {
    # generate pseudo-distrubance 'ub'
    if (method == "residual") {
      # residual bootstrapping assumes E(y|Xb) = X*b and u are iid (i.e.
      # homoskedasticity => Draw from the empirical distribution of u
      segment <- c(1:(t - p)) / (t - p)
      ub <- matrix(0, nrow = t, ncol = n)
      for (i in 1:nrow(ub)) {
        draw <- runif(1, 0, 1) # uniform distribution
        ub[i, ] <- u[min(which(segment >= draw)), ]
      }
    } else if (method == "wild") {
      # allows for heteroskedasticity
      fu <- matrix(1 - 2 * (runif(t - p, 0, 1) > 0.5), ncol = 1)
      ub <- matrix(0, nrow = (t - p), ncol = n)
      for (i in 1:n) {
        ub[, i] <- u[, i] * fu # flip sign of randomly selected 50%
      }
      ub <- rbind(matrix(0, nrow = p, ncol = n), ub)
    }

    # 2. generate pseudo-sample based on drawn u's
    # Important: The 'exogenous' interactions of endogenous variables have
    # to be constructed manually, since they directly depend on some of the
    # simulated/bootstrapped data.
    # 2a. Initial values for the artificial data: real data + plus
    # artificial residuals
    LAG <- NULL
    datab <- matrix(NA, nrow = t, ncol = n)
    for (i in 1:p) {
      datab[i, ] <- unname(as.matrix(data[i, ])) + ub[i, ]
      LAG <- c(datab[i, ], LAG)
      if (c_case == 0) {
        LAGplus <- LAG
      } else if (c_case == 1) {
        LAGplus <- c(1, LAG)
      } else if (c_case == 2) {
        LAGplus <- c(1, i, LAG)
      }
    }

    # 2b. construct interaction terms of initial observations
    if (statepos != 0) {
      if (p == 1) {
        LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos])
      } else if (p == 2) {
        LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n])
      } else if (p == 3) {
        LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n])
      } else if (p == 4) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n]
        )
      } else if (p == 5) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n]
        )
      } else if (p == 6) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n]
        )
      } else if (p == 7) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
          LAG[shockpos + 6 * n] * LAG[statepos + 6 * n]
        )
      } else if (p == 8) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
          LAG[shockpos + 6 * n] * LAG[statepos + 6 * n], LAG[shockpos + 7 * n] * LAG[statepos + 7 * n]
        )
      } else if (p == 9) {
        LAGplus <- c(
          LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
          LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
          LAG[shockpos + 6 * n] * LAG[statepos + 6 * n], LAG[shockpos + 7 * n] * LAG[statepos + 7 * n], LAG[shockpos + 8 * n] * LAG[statepos + 8 * n]
        )
      } else {
        print("Bootstrap only available up to p = 9. Please write code forward if more are needed.")
      }
    }

    #  3. generate artificial data from observation p+1 forward
    for (i in (p + 1):t) {
      for (nn in 1:n) { # forward-simulate endogenous variables
        datab[i, nn] <- LAGplus %*% A[, nn] + ub[i, nn]
      }
      LAG <- c(datab[i, ], LAG[c(1:((p - 1) * n))]) # update the LAG matrix
      if (c_case == 0) {
        LAGplus <- LAG
      } else if (c_case == 1) {
        LAGplus <- c(1, LAG)
      } else if (c_case == 2) {
        LAGplus <- c(1, i, LAG)
      }
      if (statepos != 0) { # construct the interaction variables
        if (p == 1) {
          LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos])
        } else if (p == 2) {
          LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n])
        } else if (p == 3) {
          LAGplus <- c(LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n])
        } else if (p == 4) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n]
          )
        } else if (p == 5) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n]
          )
        } else if (p == 6) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n]
          )
        } else if (p == 7) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
            LAG[shockpos + 6 * n] * LAG[statepos + 6 * n]
          )
        } else if (p == 8) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
            LAG[shockpos + 6 * n] * LAG[statepos + 6 * n], LAG[shockpos + 7 * n] * LAG[statepos + 7 * n]
          )
        } else if (p == 9) {
          LAGplus <- c(
            LAGplus, LAG[shockpos] * LAG[statepos], LAG[shockpos + n] * LAG[statepos + n], LAG[shockpos + 2 * n] * LAG[statepos + 2 * n],
            LAG[shockpos + 3 * n] * LAG[statepos + 3 * n], LAG[shockpos + 4 * n] * LAG[statepos + 4 * n], LAG[shockpos + 5 * n] * LAG[statepos + 5 * n],
            LAG[shockpos + 6 * n] * LAG[statepos + 6 * n], LAG[shockpos + 7 * n] * LAG[statepos + 7 * n], LAG[shockpos + 8 * n] * LAG[statepos + 8 * n]
          )
        }
      }
    }

    acceptdraw <- 1
    test <<- datab
    if ((sum(is.na(datab)) + sum(datab %in% Inf, na.rm = TRUE) + sum(datab %in% -Inf, na.rm = TRUE)) > 0) {
      acceptdraw <- 0
    } else if (max(abs(datab - colMeans(datab))) > 10^5) {
      acceptdraw <- 0
    } else if (qr(t(datab) %*% datab)$rank == 1) {
      acceptdraw <- 0
    } else {
      # 4. Estimate VAR, repeat if leads to errors (can happen due to
      # explosiveness of artificial data)
      if (statepos == 0) {
        exdatab <- NULL
      } else { # % if IVAR, fill in the exogenous
        exdatab <- as.matrix(datab[, shockpos] * datab[, statepos], ncol = 1)
      }
      ivar_runs <- withTimeout(
        {
          IVARb <- estimateVAR(datab, p, c_case, exdatab)
        },
        timeout = 1,
        onTimeout = "silent"
      )
      if (!is.null(ivar_runs)) { # if it failed to run within a second, skip iterataion
        IVARb <- estimateVAR(datab, p, c_case, exdatab)
      } else {
        acceptdraw <- 0
      }
    }

    # 5. Generalized impulse response functions
    if (acceptdraw == 1) {
      if (statepos == 0) {
        IVARb$GIRF1 <- simGIRF(IVARb, state, h, rep(1, length(state$s)), shocksize, fast, FALSE)
      } else { # % if IVAR, split histories in 2 and save both
        sb <- as.numeric(datab[, statepos] <= state$absval)
        sb <- sb[(p + 1):length(sb)]
        # IVARb <<- IVARb
        if (sum(sb) > 1 & sum(sb) < (length(sb) - 1)) {
          girf1_runs <- withTimeout(
            {
              IVARb$GIRF1 <- simGIRF(IVARb, state, h, !sb, shocksize, fast, FALSE)
            },
            timeout = 1,
            onTimeout = "silent"
          )
          girf2_runs <- withTimeout(
            {
              IVARb$GIRF2 <- simGIRF(IVARb, state, h, sb, shocksize, fast, FALSE)
            },
            timeout = 1,
            onTimeout = "silent"
          )
          if ((!is.null(girf1_runs)) & (!is.null(girf2_runs))) {
            IVARb$GIRF1 <- simGIRF(IVARb, state, h, !sb, shocksize, fast, FALSE)
            IVARb$GIRF2 <- simGIRF(IVARb, state, h, sb, shocksize, fast, FALSE)
            GIRF1b[, , bb] <- IVARb$GIRF1$GIRF
            GIRF2b[, , bb] <- IVARb$GIRF2$GIRF
            if ((IVARb$GIRF1$nexplos == IVARb$GIRF1$nhist) | (IVARb$GIRF2$nexplos == IVARb$GIRF2$nhist)) {
              acceptdraw <- 0
            }
          } else {
            acceptdraw <- 0
          }
        } else {
          acceptdraw <- 0
        }
      }
    }


    # 6. Decide whether or not to keep
    if (acceptdraw == 1) { # save and print progress bar
      setTxtProgressBar(pb, bb)
      bb <- bb + 1
    } else {
      nexplos <- nexplos + 1
    }
    nattempts <- nattempts + 1
  }

  # number of interrupted bootstrap iterations
  if ((nexplos / nattempts) > 0.5) {
    cat("\n Warning: ", round(100 * (nexplos / nattempts)), "% of bootstrap iterations were discarded. Possible fixes:")
    cat("\n - Start the sample at another date.")
    cat("\n - Make sure number of disregarded histories is close to zero in main GIRF function.")
    cat("\n - Choose a different bootstrap method.")
  }
  # retrieve intervals
  GIRFbands <- list()
  GIRFbands$GIRF1$GIRFbands <- array(NA, dim = c(h, n, 2))
  GIRFbands$GIRF2$GIRFbands <- array(NA, dim = c(h, n, 2))
  for (i in 1:n) {
    GIRFbands$GIRF1$GIRFbands[, i, ] <- t(apply(GIRF1b[, i, ], 1, FUN = function(x) quantile(x, probs = c(lo / 100, up / 100), na.rm = TRUE)))
    GIRFbands$GIRF2$GIRFbands[, i, ] <- t(apply(GIRF2b[, i, ], 1, FUN = function(x) quantile(x, probs = c(lo / 100, up / 100), na.rm = TRUE)))
  }

  return(GIRFbands)
}

# END