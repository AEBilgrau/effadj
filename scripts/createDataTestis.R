
################################################################################
# Creating and cleaning dataset for the paper:                                 #
#                  Variance approximation of efficiency adjusted               #
#                         DDCq-values in qPCR experiments                      #
# Script written by:                                                           #
# Steffen Falgreen, Anders Ellern Bilgrau and Martin Boegsted                  #
# Last revision: 10th of June, 2015                                            #
################################################################################

# Read metadata for testis dataset
testisMetadata <- read.csv(file = "../data/Lymphoma/testesMetadata.csv",
                           header = TRUE, sep  = ",")

# Read file names for testis dataset and exclude Yuan et al. data
# and combine in one large datafile
files <- list.files(path = "../data/Lymphoma/", pattern = "rnu|mir",
                    full.names = TRUE)

plate.number <- 1
for (file in files) {
  this.plate               <- read.csv(file = file)
  this.plate$geneName      <- gsub("\\.", "", substr(basename(file), 1, 7))
  this.plate$cDNASynthesis <-
    substr(gsub("\\.", "", substr(basename(file), 8, 100)), 1, 1)
  this.plate$cDNAMix       <- this.plate$cDNASynthesis
  this.plate$plate.number  <- plate.number

  this.plate$Well.Type <- as.character(this.plate$Well.Type)
  this.plate$Well.Type[this.plate$Well.Type == "Unknown"]  <- "sample"
  this.plate$Well.Type[this.plate$Well.Type == "Standard"] <- "standard"

  this.plate$Quantity..copies. <-
    suppressWarnings(as.numeric(as.character(this.plate$Quantity..copies.)))
  this.plate$Quantity..copies.[this.plate$Well.Type=="sample"]  <- 0
  this.plate$Quantity..copies.[this.plate$Quantity..copies.==0] <- 1

  if (plate.number == 1) {
    testis <- this.plate
  } else {
    testis <- rbind(testis, this.plate)
  }
  plate.number <- plate.number + 1
}

# Rename column-names in testis and testisMetadata

colnames(testis) <- c("wellName", "sampleName", "wellType", "threshold",
                      "Cq", "copyNumber", "geneName", "cDNAMix",
                      "cDNASynthesis", "plate.number")
testis <- merge(testis, testisMetadata)


# Determining and writing replicate number
testis$replicateNumber <- NA
replicateId  <- paste(testis$sampleName, testis$copyNumber, testis$plate.number)
replicateIds <- unique(replicateId)

for (i in replicateIds) {
  testis$replicateNumber[replicateId == i] <-
    1:length(replicateId[replicateId == i])
}

# Additional clean-up and preparation
testis$sampleType[testis$sampleType == "1"]   <- "case"
testis$sampleType[testis$sampleType == "2"]   <- "ctrl"
testis$sampleType[testis$sampleType == "NaN"] <- "Standard"
testis$sampleType <- as.factor(testis$sampleType)

testis$geneType <- "tgt"
testis$geneType[testis$geneName %in% c("rnu6b","rnu24")] <- "ref"
testis$geneType <- as.factor(testis$geneType)
testis$geneName <- as.factor(testis$geneName)

testis$threshold <- as.numeric(testis$threshold)
testis$Cq        <- suppressWarnings(as.numeric(as.character(testis$Cq)))


# Remove missing Cq-values and take only cDNA synthesis 4
# Round threshold to equal size -- should be done in MX-Pro
testis <- testis[!is.na(testis$Cq), ]
testis <- testis[testis$cDNAMix == 4 & testis$sampleName != "H420", ]
testis$threshold <- round(testis$threshold, 1)

testis$l2con <- -1*log2(testis$copyNumber)
testis$l2con <- round(testis$l2con)

# Making the data.frame into a data.qPCR class object
names(testis) <- gsub("replicateNumber", "replicate", names(testis))
class(testis) <- c("data.qPCR", "data.frame")
attr(testis, "std.curve") <- TRUE
