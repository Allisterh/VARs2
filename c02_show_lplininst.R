# c02_show_lplininst.R ----

## Settings ----
p <- 6  # easily overfitted
h <- 48
c_case <- 1
vars   <- c('GDPgap', 'Unemp', 'CoreCPIGr12', 'FFR')
exvars <- NULL
n      <- length(vars)
ident  <- 'proxy'
proxyvar  <- c('mpsprMA')
shockpos  <- 4 # position of the shock in the Cholesky ordering
shocksize <- 0 # 0 = one standard deviation, all else: absolute values
state <- list(nonlinear = 'no')
alpha <- 90 # confidence level

## Load data, generate lags, subset relevant subsample, etc. ----
load(file = 'data/data_m_ready.RData')
source('_tbx/supportfct/subset_data.R', echo = TRUE)

## Estimate matrices A, Omega, S, dynamic multipliers ----
LP <- estimateLPz(data, z, h, p, c_case, exdata, alpha, NWSE = FALSE)

## Plot impulse responses ----
plotirf1(LP$gamma, LP$gammabands, printvars)

# END
