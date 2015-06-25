
################################################################################
# Master script for the paper:                                                 #
#                                                                              #
#               Variance approximation of efficiency adjusted                  #
#                      DDCq-values in qPCR experiments                         #
#                                                                              #
#                           Bilgrau et al. (2015)                              #
#                                                                              #
# Script written by:                                                           #
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted               #
# Last revision: 10th of June, 2015                                            #
################################################################################

# Initalizing and loading required packages
if (!grepl("knitr", getwd())) setwd("./knitr")
rm(list = ls())
set.seed(987654321)
n.boots <- 1000 # Bootstrap samples in data analysis
n.sims  <- 1000 # Simulation replications

recompute <- FALSE  # Recompute heavy computations if TRUE
start.t <- proc.time()

parallel <- TRUE
n.cpus <- 4

# install.packages(c("lattice", "Hmisc", "lme4", "epiR"))
library("lattice") # For plots
library("latticeExtra")
library("lme4")    # For mixed effects models > 1.1-8
library("nlme")
library("Hmisc")   # For LaTeX tables
library("epiR")    # For performance measures and CI hereof
library("snowfall") # For parallel computing
library("GMCM")    # For multivariate simulations

# Saved R binary output
save.file <- "../output/saved.RData"
if (file.exists(save.file) && !recompute) {
  loaded <- load(save.file)
}

# 1. Load auxillary functions
source("../scripts/functions.R")

# 2. Create data sets
source("../scripts/createDataCIC.R")
source("../scripts/createDataTestis.R")
source("../scripts/createDataYuan.R")

# 3. Analyse data sets and create output
source("../scripts/createOutputCIC.R")
source("../scripts/createOutputTestis.R")
source("../scripts/createOutputYuan.R")

# 4. Run simulation experiment
source("../scripts/simulation.R")

end.t <- proc.time()
print(end.t - start.t)
