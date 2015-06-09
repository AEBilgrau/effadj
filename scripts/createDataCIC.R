
################################################################################
# Creating and cleaning dataset for the paper:                                 #
#                  Variance approximation of efficiency adjusted               #
#                         DDCq-values in qPCR experiments                      #
#              by Boegsted, Petersen, Falgreen, Bilgrau et al. (2010)          #
# Script written by:                                                           #
# Steffen Falgreen, Anders Ellern Bilgrau and Martin Boegsted                  #
# Last revision: 8th of Aug, 2012                                              #
################################################################################

#
# Read data
#

cic <- read.table(file = "../data/Testis/cic.qpcr.data.txt", header = TRUE,
                  sep  = "\t", stringsAsFactors = FALSE, dec = ",")

# Cleaning, preperations and corrections
names(cic)   <- gsub("cellLine", "sampleName", names(cic))
cic$Cq       <- suppressWarnings(as.numeric(gsub(",", ".", cic$Cq)))
cic$geneType <-
  as.factor(gsub("reference", "ref", gsub("target", "tgt", cic$geneType)))

cic$sampleType[cic$wellType == "standard"] <- "Standard"
cic$sampleType <-
  as.factor(gsub("low", "ctrl", gsub("high", "case", cic$sampleType)))

cic$sampleName <- as.factor(cic$sampleName)

cic <- cic[, !grepl("wellType", names(cic))]  # wellType is now redundant

# Remove missing Cq-values
# Round l2con values, 5.99 => 6

cic       <- cic[!is.na(cic$Cq), ]
cic$l2con <- round(cic$l2con)  # Should we round?

# Make data balanced by removing references with no paired target
cic <- subset(cic, !(cic$sampleName == "KAS-6-1" & cic$geneName == "GAPDH"))
cic <- subset(cic, !(cic$sampleName == "U-266" & cic$geneName == "GAPDH"))

cic <- subset(cic, !(cic$sampleName == "KAS-6-1" & cic$geneName == "ACTB"))
cic <- subset(cic, !(cic$sampleName == "U-266" & cic$geneName == "ACTB"))

cic <- subset(cic, !(cic$sampleName == "KAS-6-1" & cic$geneName == "MMSET"))
cic <- subset(cic, !(cic$sampleName == "U-266" & cic$geneName == "MMSET"))

#
# Converting cic into an object of class "data.qPCR";
# e.g. see the description of the "data.qPCR"-class in "Functions.R"
#

class(cic)             <- c("data.qPCR", "data.frame")
attr(cic, "std.curve") <- TRUE
