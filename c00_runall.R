rm(list = ls())

setwd("~/Dropbox/Teaching/macrometrics/codes/")
source("_tbx/supportfct/compileRfct.R")
# 1. Preliminary data stuff
source("c01_prepare_data_m.R")
source("c01_prepare_data_q.R")


# 2. Linear local projection
source("c02_show_lplinchol.R")
source("c02_show_lplininst.R")


# 2. Nonlinear local projection
source("c02_show_lpnonlin.R")
source("c02_show_lpnonlininst.R")


# 3. Linear VAR: the ancient case (Cholesky)
source("c03_show_varchol.R")


# 4. Linear VAR: sign restrictions
source("c04_show_varsign.R")


# 5. Linear VAR: IV identification
source("c05_show_proxyvar.R") # Gertler & Karadi


# 6. Smooth-transition VAR: //// only available in Matlab
# source('c06_show_stvar.R')
# source('c06_show_stvarFiscal.R')


# 7. IVAR: recession v. expansion
source("c07_show_ivar.R")
