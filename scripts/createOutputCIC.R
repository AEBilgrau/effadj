
################################################################################
# Analysis of CIC data                                                         #
# Written by:                                                                  #
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted               #
# Last revision: 11th of June, 2015                                            #
################################################################################

#
# Plot raw data
#

# Trellis plot of all data
jpeg("../output/Figure1.cic.jpg", width = 1.5*7, height = 7, units = "in",
     res = 300)
  trellis.par.set(strip.background=list(col="lightgrey"))

  fig2a <- dotplot(sampleName ~ Cq | geneType:sampleType,
                   data = subset(cic, sampleType != "Standard"),
                   main = "", col = "steelblue",
                   xlab = expression(C[q]),
                   key = list(text=list(title="A"),
                              corner = c(-0.1,1.1),
                              cex = 1.5, font = "bold"))

  fig2b <- xyplot(Cq ~ l2con | geneType:geneName, pch = 16,
                  data = subset(cic, sampleType == "Standard"),
                  xlab = as.expression(bquote(-log[2]*N["0,i,j,k"])),
                  ylab = expression(C[q]),
                  main = "", col = "steelblue",
                  key = list(text=list(title="B"),
                             corner = c(-0.1,1.1),
                             cex = 1.5, font = "bold"))

  print(fig2a, position=c(0, 0, 0.5, 1), more = TRUE)
  print(fig2b, position=c(0.5, 0, 1, 1))
dev.off()


#
# Analysis
#

#
# We wish to test:
#   MGST1 vs. GAPDH,    MGST1 vs. ACTB,      MGST1 vs. both
#   MMSET vs. GAPDH,    MMSET vs. ACTB,      MMSET vs. both
#

n.boots.cic <- n.boots
grps.list <- list(c("MGST1", "GAPDH"),
                  c("MGST1", "ACTB"),
                  c("MGST1", "ACTB", "GAPDH"),
                  c("MMSET", "GAPDH"),
                  c("MMSET", "ACTB"),
                  c("MMSET", "ACTB", "GAPDH"))

if (!exists("cic.boot") || recompute) {
  cic.boot <- list()
  for (i in 1:length(grps.list)) {
    # Subset data
    cic.tmp <- subset(cic, geneName %in% grps.list[[i]])

    # Compute bootstrap estimate
    cic.boot[[i]] <-  bootstrapEstimate(cic.tmp, n.boots = n.boots.cic)
  }
  resave(cic.boot, file = save.file)
}

# Combine results
toTeX <- NULL
for (i in 1:length(grps.list)) {
  # Subset data
  cic.tmp <- subset(cic, geneName %in% grps.list[[i]])

  results <- rbind(
    "t-test" = DDCq.test(cic.tmp, method = "N"),
    "LMEM"   = DDCq.test(cic.tmp, method = "LMM", eff.cor = F),
    "EC"     = DDCq.test(cic.tmp, method = "LMM", eff.cor = T, var.adj = F),
    "ECVA"   = DDCq.test(cic.tmp, method = "LMM", eff.cor = T, var.adj = T),
    "Bootstrap" = as.numeric(cic.boot[[i]]))

  toTeX <- rbind(toTeX, results)
}

#
# Writing LaTeX table
#

toTeX      <- signif(toTeX, 4)
toTeX[, 5] <- sn(toTeX[, 5])
colnames(toTeX) <- gsub("Pr(>|t|)", "$p$-value", colnames(toTeX), fixed = TRUE)
colnames(toTeX) <- gsub("t ", "$t$-", colnames(toTeX), fixed = TRUE)
rownames(toTeX) <- gsub("EC", "Eff. Corr.", rownames(toTeX))
rownames(toTeX) <- gsub("VA", " \\\\& Var. Adj.", rownames(toTeX))

grps <-
  sapply(grps.list, function(x) ifelse(length(x)==3,
                                       paste(x[1], "vs", x[2], "+", x[3]),
                                       paste(x[1], "vs", x[2])))

caption.txt <- "{\\bf CIC data: method comparison for estimating the
  $\\Delta\\Delta\\textrm{C}_q$-value}. $t$-test shows results from a simple
  $t$-test using only undiluted data. LMEM signifies the regular
  $\\Delta\\Delta\\textrm{C}_q$ method using a linaer mixed effects model
  without using dilution data and thus without efficiency correction.
  Eff.\\ corr.\\ denotes use of the plugin-estimator.
  Var.\\ adj.\\ denotes that the efficiency correction was variance adjusted.
  Bootstrap shows the mean and standard deviation of %d bootstrap samples."
w <- latex(toTeX,
           file    = "../output/Table1.tex",
           title   = "",
           label   = "table:cic",
           caption = sprintf(caption.txt, n.boots.cic),
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")
