
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
# Last revision: 16th of Feb., 2015                                            #
################################################################################


# Initalizing and loading required packages
rm(list = ls())
set.seed(300)

recompute <- FALSE  # Recompute heavy computations if TRUE

library("lattice")  # For plots
library("lme4")     # For mixed effects models
library("Hmisc")    # For LaTeX tables
library("ROCR")     # For ROC curves


# 1. Load auxillary functions
source("Scripts/functions.R")

# 2. Run simulation script
source("Scripts/simulation.R")

# 3. Create data sets
source("Scripts/CreateDataCIC.R", local = new.env())
source("Scripts/CreateDataTestis.R", local = new.env())

# 4. Analyse data sets and create output
source("Scripts/CreateOutputCIC.R", local = new.env())
source("Scripts/CreateOutputTestis.R", local = new.env())

# Store session information
sink(file = "Output/sessionInfo.txt")
  sessionInfo()
sink()


