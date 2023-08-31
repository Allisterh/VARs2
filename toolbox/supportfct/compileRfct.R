# ************************************************************************
# Compile R functions ----
# ************************************************************************

## Load packages ----
if (!require("gdata")) install.packages("gdata")
library(gdata)
if (!require("mFilter")) install.packages("mFilter")
library(mFilter)
if (!require("sandwich")) install.packages("sandwich")
library(sandwich)
if (!require("R.utils")) install.packages("R.utils")
library(R.utils)
if (!require("lpirfs")) install.packages("lpirfs")
library(lpirfs)

## Load own functions ----
source("toolbox/lp_tbx/estimateLP.R")
source("toolbox/lp_tbx/estimateLPz.R")
source("toolbox/lp_tbx/estimateLPnonlin.R")
source("toolbox/lp_tbx/estimateLPznonlin.R")
source("toolbox/var_tbx/estimateVAR.R")
source("toolbox/var_tbx/dyn_multipliers.R")
source("toolbox/var_tbx/bootstrapVAR.R")
source("toolbox/var_tbx/testIRFsign.R")
source("toolbox/ivar_tbx/simGIRF.R")
source("toolbox/ivar_tbx/bootstrapIVAR.R")
source("toolbox/supportfct/plotirf1.R")
source("toolbox/supportfct/plotirf2.R")
source("toolbox/supportfct/plotirf3.R")

## Other settings ----
par(mar = c(4.6, 4.1, 2.1, 2.1))

# END