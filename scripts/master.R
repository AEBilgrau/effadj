
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
# setwd("./knitr")
rm(list = ls())
set.seed(987654321)

recompute <- FALSE  # Recompute heavy computations if TRUE

# library("Bmisc")   # Resave function and more http://github.com/AEBilgrau/Bmisc
library("lattice") # For plots
library("lme4")    # For mixed effects models
library("Hmisc")   # For LaTeX tables
library("ROCR")    # For ROC curves

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

# 2. Run simulation script
source("../scripts/simulation.R")

# 3. Create data sets
source("../scripts/createDataCIC.R")
source("../scripts/createDataTestis.R")

# 4. Analyse data sets and create output
source("../scripts/createOutputCIC.R", local = new.env())
source("../scripts/createOutputTestis.R", local = new.env())

