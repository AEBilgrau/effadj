
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

library("Bmisc")    # Resave and more. http://github.com/AEBilgrau/Bmisc
library("lattice")  # For plots
library("lme4")     # For mixed effects models
library("Hmisc")    # For LaTeX tables
library("ROCR")     # For ROC curves

# Saved R binary output
save.file <- "output/saved.RData"
if (file.exists(save.file) | recompute) {
  load(save.file)
}

# Colors
jet.colors <-
  colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
                     "#7FFF7F", "yellow", "#FF7F00", "red"))


# 1. Load auxillary functions
source("Scripts/functions.R")

# 2. Run simulation script
source("Scripts/simulation.R")

# 3. Create data sets
source("Scripts/createDataCIC.R", local = new.env())
source("Scripts/createDataTestis.R", local = new.env())

# 4. Analyse data sets and create output
source("Scripts/createOutputCIC.R", local = new.env())
source("Scripts/createOutputTestis.R", local = new.env())


