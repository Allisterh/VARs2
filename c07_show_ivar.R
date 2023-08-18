# c07_show_ivar.R ----

## Settings ----
p <-  4
h <- 25
c_case <- 1
vars   <- c('GDPgap','Unemp','CoreCPIGr4','FFRshadow')
exvars <- NULL  # has to be empty in IVAR - will be constructed in subset_data.R
n      <- length(vars)
ident  <- 'chol'
shockpos  <- 4 # position of the shock in the Cholesky ordering
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear  = 'yes',
              logistic   = 'no',
              interacted = 'yes',
              shockvar   = 'FFRshadow',  # both 'shockvar' and 'shockpos' need to be part of 'vars'
              statevar   = 'GDPgap',
              cq = 12)
nboot <- 1000
alpha <-   90 # confidence level

## Load data, generate lags, subset relevant subsample, etc. ----
load('data/data_q_ready.RData') # quarterly!
source('_tbx/supportfct/subset_data.R', echo = TRUE)

## Estimate a linear VAR with interaction term in exdata ----
IVAR <- estimateVAR(data, p, c_case, exdata)

## Split sample into histories of two groups and simulate generalized impulse responses
IVAR$GIRF1 <- simGIRF(IVAR, state, h, !state$s[(p+1):length(state$s)], shocksize, TRUE, TRUE) # do h + a number to make sure GIRFs are stable
IVAR$GIRF2 <- simGIRF(IVAR, state, h, state$s[(p+1):length(state$s)], shocksize, TRUE, TRUE)

## Bootstrap ----
IVAR$GIRFbands <- bootstrapIVAR(IVAR, state, nboot, alpha, 'residual', TRUE)

## Figure ----
plotirf2(IVAR$GIRF1$GIRF[1:h, ], IVAR$GIRFbands$GIRF1$GIRFbands[1:h, , ], 
         IVAR$GIRF2$GIRF[1:h, ], IVAR$GIRFbands$GIRF2$GIRFbands[1:h, , ], 
         printvars, c('Boom', 'Recession'))

# END