# ************************************************************************
# Simulate GIRF ----
# ************************************************************************
simGIRF <- function(IVAR, state, h, historyindex, impulse, fast = TRUE, printprogress = TRUE) {
  # references Koop, Peseran & Potter (1996), Giovanni Pellegrino
  # OIRFhist(t,j,k,hist) : GIRFs associated to each starting history in the input vector "histories"

  ## Unpack and re-define settings
  c_case <- IVAR$c_case
  n <- IVAR$n
  p <- IVAR$p
  A <- IVAR$A
  S <- IVAR$S
  uu <- IVAR$u
  historyindex <- c(FALSE, historyindex[1:(length(historyindex) - 1)]) # lag by 1 period
  histories <- cbind(IVAR$X, IVAR$Xex)
  if (c_case == 1) {
    histories <- cbind(
      matrix(1, nrow = nrow(histories), ncol = 1),
      histories
    )
  } else if (c_case == 2) {
    histories <- cbind(
      matrix(1, nrow = nrow(histories), ncol = 1),
      matrix(seq(1, nrow(histories), 1), nrow = nrow(histories), ncol = 1),
      histories
    )
  }
  histories <- histories[historyindex == 1, ]
  nhist <- nrow(histories)
  npx <- ncol(histories) # % number of coefficients for each equation: = c_case + n*p (endogenous lags) + p (exogenous lags)

  shockpos <- state$shockpos
  statepos <- state$statepos

  iter <- 500 # number of repetitions simulated for each initial condition
  if (fast) {
    iter <- 2 # no extraction of future shocks, 2 draws sufficient
  }
  mode <- 0 # 0 = draw from empirical distribution of residuals

  nexplos <- 0 # to count for the number of explosive histories
  nrep <- rep(0, nhist)
  IRFhist <- array(0, dim = c(h, n, nhist))
  if (printprogress) {
    print("Progress...")
    pb <- txtProgressBar(min = 0, max = nhist, initial = 0, style = 3)
  }

  ## Loop over all histories defined in input
  for (ii in 1:nhist) {
    if (printprogress) {
      setTxtProgressBar(pb, ii)
    }

    irf <- array(0, dim = c(h, n, iter)) # lower case irf: h x n x draws (500)

    ic <- 1
    while (ic <= iter) {
      # preparations, containers
      u <- matrix(0, nrow = h, ncol = n)
      eps <- matrix(0, nrow = 1, ncol = n)
      ushocked <- matrix(0, nrow = 1, ncol = n)
      X <- histories[ii, ] # predefined matrix for the starting values for the RHS variables
      W <- histories[ii, ]
      Y <- matrix(0, nrow = h, ncol = n) # variables including shock
      Z <- matrix(0, nrow = h, ncol = n) # variable excluding shock

      # settings
      if (!fast) {
        if (mode == 0) {
          u <- uu[ceiling(nrow(uu) * runif(h)), ] # for each h of horizon and variable, draw shock from empirical distribution
        } else {
          u <- matrix(rnorm(h * n), nrow = n) %*% t(S)
        }
      }
      eps[1, ] <- solve(S) %*% t(t(u[1, ])) # orthogonalized shock

      # add an additional shock to the structural eps vector
      if (impulse == 0) {
        shock <- 1
        eps[1, shockpos] <- eps[1, shockpos] + shock
      } else {
        shock <- impulse / S[shockpos, shockpos] # initial shock (in terms of st.dev.)
        eps[1, shockpos] <- eps[1, shockpos] + shock
      }
      ushocked[1, ] <- S %*% t(t(eps[1, ])) # re-transform to reduced-form residual

      # period 0
      Y[1, ] <- X %*% A + ushocked[1, ]
      Z[1, ] <- W %*% A + u[1, ] # only usual residuals
      irf[1, , ic] <- Y[1, ] - Z[1, ] # response in the period 0

      # period 1 to h-1
      for (tt in 1:(h - 1)) {
        NewX <- rep(0, npx)
        NewW <- rep(0, npx)
        if (c_case == 2) { # if c_case = 2, then RHS includes the time
          time <- histories[ii, 2]
          NewX[1:2] <- c(1, time + tt)
          NewW[1:2] <- c(1, time + tt)
        } else if (c_case == 1) {
          NewX[1] <- 1
          NewW[1] <- 1
        }
        # the amount by which X (and W) are 'moved' to the left every iteration depends on the number of lags (p)
        if (p == 1) {
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 2) {
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 3) {
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 4) {
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 5) {
          NewX[(1 + c_case + 4 * n):(c_case + 5 * n)] <- X[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewW[(1 + c_case + 4 * n):(c_case + 5 * n)] <- W[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 6) {
          NewX[(1 + c_case + 5 * n):(c_case + 6 * n)] <- X[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewW[(1 + c_case + 5 * n):(c_case + 6 * n)] <- W[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewX[(1 + c_case + 4 * n):(c_case + 5 * n)] <- X[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewW[(1 + c_case + 4 * n):(c_case + 5 * n)] <- W[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 7) {
          NewX[(1 + c_case + 6 * n):(c_case + 7 * n)] <- X[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewW[(1 + c_case + 6 * n):(c_case + 7 * n)] <- W[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewX[(1 + c_case + 5 * n):(c_case + 6 * n)] <- X[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewW[(1 + c_case + 5 * n):(c_case + 6 * n)] <- W[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewX[(1 + c_case + 4 * n):(c_case + 5 * n)] <- X[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewW[(1 + c_case + 4 * n):(c_case + 5 * n)] <- W[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 8) {
          NewX[(1 + c_case + 7 * n):(c_case + 8 * n)] <- X[(1 + c_case + 6 * n):(c_case + 7 * n)]
          NewW[(1 + c_case + 7 * n):(c_case + 8 * n)] <- W[(1 + c_case + 6 * n):(c_case + 7 * n)]
          NewX[(1 + c_case + 6 * n):(c_case + 7 * n)] <- X[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewW[(1 + c_case + 6 * n):(c_case + 7 * n)] <- W[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewX[(1 + c_case + 5 * n):(c_case + 6 * n)] <- X[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewW[(1 + c_case + 5 * n):(c_case + 6 * n)] <- W[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewX[(1 + c_case + 4 * n):(c_case + 5 * n)] <- X[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewW[(1 + c_case + 4 * n):(c_case + 5 * n)] <- W[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else if (p == 9) {
          NewX[(1 + c_case + 8 * n):(c_case + 9 * n)] <- X[(1 + c_case + 7 * n):(c_case + 8 * n)]
          NewW[(1 + c_case + 8 * n):(c_case + 9 * n)] <- W[(1 + c_case + 7 * n):(c_case + 8 * n)]
          NewX[(1 + c_case + 7 * n):(c_case + 8 * n)] <- X[(1 + c_case + 6 * n):(c_case + 7 * n)]
          NewW[(1 + c_case + 7 * n):(c_case + 8 * n)] <- W[(1 + c_case + 6 * n):(c_case + 7 * n)]
          NewX[(1 + c_case + 6 * n):(c_case + 7 * n)] <- X[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewW[(1 + c_case + 6 * n):(c_case + 7 * n)] <- W[(1 + c_case + 5 * n):(c_case + 6 * n)]
          NewX[(1 + c_case + 5 * n):(c_case + 6 * n)] <- X[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewW[(1 + c_case + 5 * n):(c_case + 6 * n)] <- W[(1 + c_case + 4 * n):(c_case + 5 * n)]
          NewX[(1 + c_case + 4 * n):(c_case + 5 * n)] <- X[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewW[(1 + c_case + 4 * n):(c_case + 5 * n)] <- W[(1 + c_case + 3 * n):(c_case + 4 * n)]
          NewX[(1 + c_case + 3 * n):(c_case + 4 * n)] <- X[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewW[(1 + c_case + 3 * n):(c_case + 4 * n)] <- W[(1 + c_case + 2 * n):(c_case + 3 * n)]
          NewX[(1 + c_case + 2 * n):(c_case + 3 * n)] <- X[(1 + c_case + n):(c_case + 2 * n)]
          NewW[(1 + c_case + 2 * n):(c_case + 3 * n)] <- W[(1 + c_case + n):(c_case + 2 * n)]
          NewX[(1 + c_case + n):(c_case + 2 * n)] <- X[(1 + c_case):(c_case + n)]
          NewW[(1 + c_case + n):(c_case + 2 * n)] <- W[(1 + c_case):(c_case + n)]
          NewX[(1 + c_case):(c_case + n)] <- Y[tt, ]
          NewW[(1 + c_case):(c_case + n)] <- Z[tt, ]
        } else {
          print("GIRFs only available up to p = 9. Please write code forward if more are needed.")
        }

        # Based on the re-computed endogenous variables, we re-compute
        # the 'exogenous' interactions
        if (statepos != 0) {
          if (p == 1) {
            ordint <- c_case + shockpos
            ordint2 <- c_case + statepos
            NewX[npx] <- NewX[ordint2] * NewX[ordint]
            NewW[npx] <- NewW[ordint2] * NewW[ordint]
          } else if (p == 2) {
            ordint <- c(c_case + shockpos, c_case + n + shockpos)
            ordint2 <- c(c_case + statepos, c_case + n + statepos)
            NewX[npx - 1] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewW[npx - 1] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx] <- NewW[ordint2[2]] * NewW[ordint[2]]
          } else if (p == 3) {
            ordint <- c(c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos)
            ordint2 <- c(c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos)
            NewX[npx - 2] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 1] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewW[npx - 2] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 1] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx] <- NewW[ordint2[3]] * NewW[ordint[3]]
          } else if (p == 4) {
            ordint <- c(c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos)
            ordint2 <- c(c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos)
            NewX[npx - 3] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 2] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 1] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewW[npx - 3] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 2] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 1] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx] <- NewW[ordint2[4]] * NewW[ordint[4]]
          } else if (p == 5) {
            ordint <- c(
              c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos,
              c_case + 4 * n + shockpos
            )
            ordint2 <- c(
              c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos,
              c_case + 4 * n + statepos
            )
            NewX[npx - 4] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 3] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 2] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx - 1] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewX[npx] <- NewX[ordint2[5]] * NewX[ordint[5]]
            NewW[npx - 4] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 3] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 2] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx - 1] <- NewW[ordint2[4]] * NewW[ordint[4]]
            NewW[npx] <- NewW[ordint2[5]] * NewW[ordint[5]]
          } else if (p == 6) {
            ordint <- c(
              c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos,
              c_case + 4 * n + shockpos, c_case + 5 * n + shockpos
            )
            ordint2 <- c(
              c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos,
              c_case + 4 * n + statepos, c_case + 5 * n + statepos
            )
            NewX[npx - 5] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 4] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 3] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx - 2] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewX[npx - 1] <- NewX[ordint2[5]] * NewX[ordint[5]]
            NewX[npx] <- NewX[ordint2[6]] * NewX[ordint[6]]
            NewW[npx - 5] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 4] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 3] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx - 2] <- NewW[ordint2[4]] * NewW[ordint[4]]
            NewW[npx - 1] <- NewW[ordint2[5]] * NewW[ordint[5]]
            NewW[npx] <- NewW[ordint2[6]] * NewW[ordint[6]]
          } else if (p == 7) {
            ordint <- c(
              c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos,
              c_case + 4 * n + shockpos, c_case + 5 * n + shockpos, c_case + 6 * n + shockpos
            )
            ordint2 <- c(
              c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos,
              c_case + 4 * n + statepos, c_case + 5 * n + statepos, c_case + 6 * n + statepos
            )
            NewX[npx - 6] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 5] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 4] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx - 3] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewX[npx - 2] <- NewX[ordint2[5]] * NewX[ordint[5]]
            NewX[npx - 1] <- NewX[ordint2[6]] * NewX[ordint[6]]
            NewX[npx] <- NewX[ordint2[7]] * NewX[ordint[7]]
            NewW[npx - 6] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 5] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 4] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx - 3] <- NewW[ordint2[4]] * NewW[ordint[4]]
            NewW[npx - 2] <- NewW[ordint2[5]] * NewW[ordint[5]]
            NewW[npx - 1] <- NewW[ordint2[6]] * NewW[ordint[6]]
            NewW[npx] <- NewW[ordint2[7]] * NewW[ordint[7]]
          } else if (p == 8) {
            ordint <- c(
              c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos,
              c_case + 4 * n + shockpos, c_case + 5 * n + shockpos, c_case + 6 * n + shockpos,
              c_case + 7 * n + shockpos
            )
            ordint2 <- c(
              c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos,
              c_case + 4 * n + statepos, c_case + 5 * n + statepos, c_case + 6 * n + statepos,
              c_case + y * n + statepos
            )
            NewX[npx - 7] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 6] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 5] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx - 4] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewX[npx - 3] <- NewX[ordint2[5]] * NewX[ordint[5]]
            NewX[npx - 2] <- NewX[ordint2[6]] * NewX[ordint[6]]
            NewX[npx - 1] <- NewX[ordint2[7]] * NewX[ordint[7]]
            NewX[npx] <- NewX[ordint2[8]] * NewX[ordint[8]]
            NewW[npx - 7] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 6] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 5] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx - 4] <- NewW[ordint2[4]] * NewW[ordint[4]]
            NewW[npx - 3] <- NewW[ordint2[5]] * NewW[ordint[5]]
            NewW[npx - 2] <- NewW[ordint2[6]] * NewW[ordint[6]]
            NewW[npx - 1] <- NewW[ordint2[7]] * NewW[ordint[7]]
            NewW[npx] <- NewW[ordint2[8]] * NewW[ordint[8]]
          } else if (p == 9) {
            ordint <- c(
              c_case + shockpos, c_case + n + shockpos, c_case + 2 * n + shockpos, c_case + 3 * n + shockpos,
              c_case + 4 * n + shockpos, c_case + 5 * n + shockpos, c_case + 6 * n + shockpos,
              c_case + 7 * n + shockpos, c_case + 8 * n + shockpos
            )
            ordint2 <- c(
              c_case + statepos, c_case + n + statepos, c_case + 2 * n + statepos, c_case + 3 * n + statepos,
              c_case + 4 * n + statepos, c_case + 5 * n + statepos, c_case + 6 * n + statepos,
              c_case + y * n + statepos, c_case + 8 * n + statepos
            )
            NewX[npx - 8] <- NewX[ordint2[1]] * NewX[ordint[1]]
            NewX[npx - 7] <- NewX[ordint2[2]] * NewX[ordint[2]]
            NewX[npx - 6] <- NewX[ordint2[3]] * NewX[ordint[3]]
            NewX[npx - 5] <- NewX[ordint2[4]] * NewX[ordint[4]]
            NewX[npx - 4] <- NewX[ordint2[5]] * NewX[ordint[5]]
            NewX[npx - 3] <- NewX[ordint2[6]] * NewX[ordint[6]]
            NewX[npx - 2] <- NewX[ordint2[7]] * NewX[ordint[7]]
            NewX[npx - 1] <- NewX[ordint2[8]] * NewX[ordint[8]]
            NewX[npx] <- NewX[ordint2[9]] * NewX[ordint[9]]
            NewW[npx - 8] <- NewW[ordint2[1]] * NewW[ordint[1]]
            NewW[npx - 7] <- NewW[ordint2[2]] * NewW[ordint[2]]
            NewW[npx - 6] <- NewW[ordint2[3]] * NewW[ordint[3]]
            NewW[npx - 5] <- NewW[ordint2[4]] * NewW[ordint[4]]
            NewW[npx - 4] <- NewW[ordint2[5]] * NewW[ordint[5]]
            NewW[npx - 3] <- NewW[ordint2[6]] * NewW[ordint[6]]
            NewW[npx - 2] <- NewW[ordint2[7]] * NewW[ordint[7]]
            NewW[npx - 1] <- NewW[ordint2[8]] * NewW[ordint[8]]
            NewW[npx] <- NewW[ordint2[9]] * NewW[ordint[9]]
          }
        }

        X <- NewX
        W <- NewW
        Y[tt + 1, ] <- X %*% A + u[tt + 1, ] # next period's endogenous variables
        Z[tt + 1, ] <- W %*% A + u[tt + 1, ]
        irf[tt + 1, , ic] <- Y[tt + 1, ] - Z[tt + 1, ] # %store the empirical response fot the period tt
      }

      # repeat the draw if residuals extracted make the GIRF explosive
      if ((sum(is.na(irf)) == 0) & (sum(is.na(S)) == 0) & (irf[h, shockpos, ic] <= S[shockpos, shockpos] * 10 * abs(shock)) & (irf[h, shockpos, ic] >= -S[shockpos, shockpos] * 10 * abs(shock))) {
        ic <- ic + 1 # % not explosive: move to next draw
      } else { #  the draw is repeated since ic remain the same...
        nrep[ii] <- nrep[ii] + 1
        if (nrep[ii] == iter / 2) { # ...but if the same extraction has been repeated for
          # much of the iter considered it means that is not only due to a particular sequence
          # of future shocks, but rather that the starting history considered is explosive in itself
          nexplos <- nexplos + 1
          irf <- array(0, dim = c(h, n, iter)) # to not count at all when computing the average
          next # to leave the while loop and consider the next history
        }
      }
    } # end loop over iterations

    IRFhist[, , ii] <- apply(irf, c(1, 2), mean) # mean across iterations of a particular history
  }

  if (printprogress) {
    cat("\n Number of GIRFs that are explosive: ", nexplos, "/", nhist, "\n")
  }

  ## Average the mean response across (non-explosive) histories
  IRFhistsum <- matrix(0, nrow = h, ncol = n)
  for (hh in 1:nhist) {
    IRFhistsum <- IRFhistsum + IRFhist[, , hh]
  }
  GIRF <- IRFhistsum / (nhist - nexplos) # % sample mean across the non-explosive IRF

  ## Output
  return(GIRF = list(
    "GIRF" = GIRF,
    "S" = S,
    "h" = h,
    "statepos" = statepos,
    "shockpos" = shockpos,
    "impulse" = impulse,
    "nexplos" = nexplos,
    "nhist" = nhist,
    "nrep" = nrep
  ))
}

# END