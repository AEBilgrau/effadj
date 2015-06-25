################################################################################
# Yuan analysis                                                                #
# Script written by:                                                           #
#   Anders Ellern Bilgrau                                                      #
################################################################################

#
# Download and extract Yuan et al 2008 data
#

yuan.destfile <- "../data/Yuan/200700169_s.zip"
if (!file.exists(yuan.destfile)) {
  download.file("http://www.wiley-vch.de/contents/jc_2446/2008/200700169_s.zip",
                destfile = yuan.destfile, method = "internal", mode = "wb")
  unzip(yuan.destfile, exdir = dirname(yuan.destfile))
}

# Read table
yuan <- read.table(file = "../data/Yuan/200700169_s8.txt",
                   col.names = c("treat","geneName","copyNumber","Cq","group"))

# Organize columns
yuan$sampleType <- factor(yuan$treat,
                          levels = c("Ala2h", "water2h"),
                          labels = c("case", "ctrl"))
yuan$geneType <- factor(as.character(yuan$geneName))
levels(yuan$geneType) <- c("tgt", "ref", "ref")
yuan$l2con <- -log2(yuan$copyNumber/25)
yuan$replicate <- 1

# Create sample names
s <- 1
for (i in unique(yuan$sampleType)) {
  for (l in unique(yuan$l2con)) {
    get <- with(yuan, sampleType == i & l2con == l)
    yuan$sampleName[get] <- rep(sprintf("s%04d", s), sum(get))
    s <- s + 1
  }
}

# Aggregate, mean over technical replicates
# yuan <- aggregate(Cq ~ ., FUN = mean, data = yuan)


