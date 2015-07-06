
################################################################################
# Analysis of CIC data                                                         #
# Written by:                                                                  #
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted               #
# Last revision: 11th of June, 2015                                            #
################################################################################

#
# Plot raw data
#

mypanel <- function(x, y, ...) {

  panel.lmlineq(x, y, adj = c(1,0), lty = 2, col.line = "darkgrey", digits = 2,
                at = 0.5, rot = TRUE, cex = 0.7, style = 2, pos = 1,
                varNames = alist(y = C[q], x = k), offset = 1)
  panel.xyplot(x, y, ...)
  aa <- lm(y ~ x)
  panel.text(2, 31.5, labels = sprintf("se = %.3f", sqrt(vcov(aa)[2,2])),
             cex = 0.7)
}

# Trellis plot of all data
library("ggplot2")
blue.colours <- c("#5C7495", "#4D6FC0", "#3E424B", "#2D406A")

cic.data <- subset(cic, l2con == 0)
levels(cic.data$sampleType) <-
  gsub("Standard", "case", levels(cic.data$sampleType))

fig1a <-
  dotplot(sampleName ~ Cq | geneType:sampleType,
          groups = cic.data$geneName,
          col = blue.colours,
          data = cic.data,
          main = "",
          xlab = expression(C[q]),
          pch = seq_along(blue.colours),
          key = list(text = list(title = "A"),
                     corner = c(-0.1,1.1),
                     cex = 1.5, font = "bold"))

cic.std <- subset(cic, sampleName == "AMO-1")
fig1b <-
  xyplot(Cq ~ l2con | as.factor(geneName),
         groups = cic.std$geneName,
         col = blue.colours,
         panel = mypanel,
         data = cic.std,
         xlab = as.expression(bquote(-log[2]*N["0,i,j,k"])),
         ylab = expression(C[q]),
         pch = seq_along(blue.colours),
         main = "",
         key = list(text = list(title = "B"),
                    corner = c(-0.1,1.1),
                    cex = 1.5, font = "bold"))

setEPS()
postscript("../output/fig1.eps", width = 1.5*7, height = 0.5*7, fonts = "serif")
  trellis.par.set(strip.background = list(col = "lightgrey"))
  print(fig1a, position = c(0, 0, 0.5, 1), more = TRUE)
  print(fig1b, position = c(0.5, 0, 1, 1))
dev.off()

rm(cic.data, cic.std)

#
# Analysis
#

#
# We wish to test:
#   MGST1 vs. GAPDH,    MGST1 vs. ACTB,      MGST1 vs. both (omitted)
#   MMSET vs. GAPDH,    MMSET vs. ACTB,      MMSET vs. both (omitted)
#

grps.list.cic <- list(c("MGST1", "GAPDH"),
                      c("MGST1", "ACTB"),
                      c("MMSET", "GAPDH"),
                      c("MMSET", "ACTB"))


if (!exists("cic.boot") || recompute) {
  message("CIC boostrap")
  cic.boot <- list()
  for (i in 1:length(grps.list.cic)) {
    # Subset data
    cic.tmp <- as.data.qPCR(subset(cic, geneName %in% grps.list.cic[[i]]))

    # Compute bootstrap estimate
    cic.boot[[i]]  <- bootstrapEstimate(cic.tmp, n.boots = n.boots,
                                        weighted = TRUE, alpha = 0.05)

    message(sprintf("i = %d", i))
  }

  resave(cic.boot, file = save.file)
}

# Combine results
sink("../output/fit_summary_cic.txt")
pdf("../output/modelchecks_cic.pdf", onefile = TRUE)
toTeX <- NULL
we <- TRUE
for (i in 1:length(grps.list.cic)) {
  # Subset data
  cic.tmp <- subset(cic, geneName %in% grps.list.cic[[i]])

  # Print fit information
  m <- paste(grps.list.cic[[i]], collapse = " vs. ")
  cat("\n\n\n\n\n===", m,"===\n")
  print(summary(fit <- qPCRfit(as.data.qPCR(cic.tmp), weighted = we)))
  print(plot(fit, col = cic.tmp$sampleName, pch = cic.tmp$sampleType),
        main = m)

  # Create results for table
  results <- rbind(
    "t-test" = DDCq.test(cic.tmp, method = "N"),
    "LMM"    = DDCq.test(cic.tmp, method = "LMM",
                         eff.cor = F, weighted = we),
    "EC"     = DDCq.test(cic.tmp, method = "LMM",
                         eff.cor = T, var.adj = F, weighted = we),
    "ECVA1"  = DDCq.test(cic.tmp, method = "LMM",
                         eff.cor = T, var.adj = T, weighted = we),
    "ECVA2"  = DDCq.test(cic.tmp, method = "LMM",
                         eff.cor = T, var.adj = T,
                         var.type = "montecarlo", weighted = we),
    "Bootstrap" = as.numeric(cic.boot[[i]])
    )

  toTeX <- rbind(toTeX, results)
}
dev.off()
sink()



#
# Writing LaTeX table
#

toTeX      <- signif(toTeX, 4)
toTeX[, 5] <- sn(toTeX[, 5])
colnames(toTeX) <- gsub("Pr(>|t|)", "$p$-value", colnames(toTeX), fixed = TRUE)
colnames(toTeX) <- gsub("t ", "$t$-", colnames(toTeX), fixed = TRUE)

rownames(toTeX) <- gsub("LMEM", "LMM", rownames(toTeX))
rownames(toTeX) <- gsub("t.", "$t$-", rownames(toTeX), fixed = TRUE)
rownames(toTeX) <- gsub("ECVA", "EC\\\\&VA", rownames(toTeX))

grps <-
  sapply(grps.list.cic, function(x) ifelse(length(x) == 3,
                                           paste(x[1], "vs", x[2], "+", x[3]),
                                           paste(x[1], "vs", x[2])))

caption.txt <- "CIC data: Method comparison for estimating the
  $\\ddcq$-value.
  EC denotes use of the plugin-estimator.
  VA denotes that the efficiency correction was variance adjusted using the
  delta method (1) or Monte Carlo integration (2).
  Bootstrap shows the mean and standard deviation of %d
  bootstrap samples using the EC estimate. The last two columns show the $95%s$
  lower and upper confidence interval limits."

toTeX <- toTeX[!grepl("LMM|t-test", rownames(toTeX)), ]

w <- latex(toTeX,
           file    = "../output/Table1.tex",
           title   = "",
           label   = "table:cic",
           caption = sprintf(caption.txt, n.boots, "\\%"),
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")
