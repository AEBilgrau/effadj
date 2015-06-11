
################################################################################
# Master script for the paper:                                                 #
#                                                                              #
#               Variance approximation of efficiency adjusted                  #
#                      DDCq-values in qPCR experiments                         #
#                                                                              #
#                           Bilgrau et al. (2011)                              #
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
n.sims  <- 100  # Simulation replications

recompute <- FALSE  # Recompute heavy computations if TRUE

# install.packages(c("lattice", "Hmisc", "lme4", "epiR"))
library("lattice") # For plots
library("lme4")    # For mixed effects models
library("Hmisc")   # For LaTeX tables
library("epiR")    # For performance measures and CI hereof

# Saved R binary output
save.file <- "../output/saved.RData"
if (file.exists(save.file) || recompute) {
  load(save.file)
}


# Colors
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                                 "#7FFF7F", "yellow", "#FF7F00", "red"))

# 1. Load auxillary functions
source("../scripts/functions.R")

# 2. Create data sets
source("../scripts/createDataCIC.R")
source("../scripts/createDataTestis.R")

# 3. Analyse data sets and create output
source("../scripts/createOutputCIC.R")
source("../scripts/createOutputTestis.R")

# 4. Run simulation experiment
source("../scripts/simulation.R")
