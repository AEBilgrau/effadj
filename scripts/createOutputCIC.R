
################################################################################
# Creates LaTeX tables                                                         #
# Script written by:                                                           #
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted               #
# Last revision: 10th of June, 2013                                            #
################################################################################

# Trellis plot of all data
jpeg("../output/Figure2.cic.jpg", width = 7*2, height = 4, units = "in",
     res = 300)
  trellis.par.set(strip.background=list(col="lightgrey"))

  fig2a <- dotplot(sampleName ~ Cq | geneType:sampleType,
                   data = cic[cic$sampleType != "Standard", ],
                   main = "", col = "steelblue",
                   key = list(text=list(title="A"),
                              corner = c(-0.1,1.1),
                              cex = 1.5, font = "bold"))

  fig2b <- dotplot(l2con ~ Cq | geneType:factor(geneName),
                   scales = list(y = list(lab = sort(unique(cic$l2con)))),
                   data = cic[cic$sampleType == "Standard", ],
                   ylab = as.expression(bquote(-log[2]*N["0,i,j,k"])),
                   main = "", col = "steelblue",
                   key = list(text=list(title="B"),
                              corner = c(-0.1,1.1),
                              cex = 1.5, font = "bold"))

  print(fig2a, position=c(0, 0, 0.5, 1), more = TRUE)
  print(fig2b, position=c(0.5, 0, 1, 1))
dev.off()


#
# Bootstrapped confidence bounds and p-values
#

new.cic <- aggregate(Cq ~ sampleName + geneType + sampleType +
                            geneType + l2con + copyNumber,
                     data = cic, FUN = mean)
new.cic <- as.data.qPCR(new.cic)  # Making new.cic a data.qPCR object

start.sample   <- 3
n.samples      <- 6
start.dilution <- 3
n.dilutions    <- 5
end.sample     <- start.sample + n.samples - 1

if (!exists("cic.boot") || recompute) {
  cic.boot <- Bootstrap.qPCR(new.cic,
                             start.sample   = start.sample,
                             n.samples      = n.samples,
                             start.dilution = start.dilution,
                             n.dilutions    = n.dilutions,
                             n.resamp       = 1000)
  resave(cic.boot, file = save.file)
}

# Plotting bootstrap results
pdf("../output/cic.pcurve.pdf")
  plot(1, type = "n",
       xlab = "Samples per group",
       xlim = c(start.sample, start.sample + n.samples - 1),
       ylab = "P-value",
       ylim = c(0, 1), log = "")
  grid()

  lines(c(start.sample, end.sample), c(0.05, 0.05), lty = 2)

  for (i in 1:length(cic.boot)) {
    lines(start.sample:end.sample, cic.boot[[i]][4, ],
          type = "b", col = i, lwd = 1.5)
    lines(start.sample:end.sample, cic.boot[[i]][5, ],
          type = "b", col = i, lwd = 1.5, lty = 2)
    lines(start.sample:end.sample, cic.boot[[i]][6, ],
          type = "b", col = i, lwd = 1.5, lty = 2)
  }
  legend("right",
         legend = paste(start.dilution + 0:(n.dilutions - 1), "step dilution"),
         lty = 1, lwd = 1.5, bty = "n", title = "", inset = 0.05,
         col = 1:length(cic.boot))
dev.off()

pdf("../output/cic.estcurve.pdf")
  plot(1, type = "n",
       xlab = "Number of samples per group",
       xlim = c(start.sample, start.sample + n.samples - 1),
       ylab = "DDCq",
       ylim = range(cic.boot) + 0.01*c(-1,1))
  grid()

  for (i in 1:length(cic.boot)) {
    lines(start.sample:end.sample, cic.boot[[i]][1,],
          type = "b", col = i, lwd = 1.5)
    lines(start.sample:end.sample, cic.boot[[i]][2,],
          type = "b", col = i, lwd = 1.5, lty = 2)
    lines(start.sample:end.sample, cic.boot[[i]][3,],
          type = "b", col = i, lwd = 1.5, lty = 2)
  }
  legend("topright",
         legend = paste(start.dilution + 0:(n.dilutions - 1), "step dilution"),
         lty = 1, lwd = 1.5, bty = "n", title = "",
         col = 1:length(cic.boot))
dev.off()


#
# Analysis
#

#
# We wish to test:
# MGST1 vs. GAPDH,  MGST1 vs. ACTB,  MGST1 vs. both
# MMSET vs. GAPDH,  MMSET vs. ACTB,  MMSET vs. both
#

grps.list <- list(c("MGST1", "GAPDH"),
                  c("MGST1", "ACTB"),
                  c("MGST1", "ACTB", "GAPDH"),
                  c("MMSET", "GAPDH"),
                  c("MMSET", "ACTB"),
                  c("MMSET", "ACTB", "GAPDH"))

toTeX <- NULL
for (i in 1:length(grps.list)) {
  genes <- grps.list[[i]]

  cic.tmp <- as.data.qPCR(cic[cic$geneName %in% genes, ])
  results <-
    rbind("t-test" = DDCq.test(cic.tmp, method = "N"),
          "LMEM"   = DDCq.test(cic.tmp, method = "LMM", eff.cor = FALSE),
          "EC"     = DDCq.test(cic.tmp, method = "LMM",
                               eff.cor = TRUE, var.adj = FALSE),
          "ECVA"   = DDCq.test(cic.tmp, method = "LMM",
                               eff.cor = TRUE, var.adj = TRUE))

  toTeX <- rbind(toTeX, results)
}

#
# Writing LaTeX table
#

toTeX      <- signif(toTeX, 4)
toTeX[, 5] <- sn(toTeX[, 5])
colnames(toTeX) <- gsub(">|t|", "$>|t|$", colnames(toTeX), fixed = TRUE)
rownames(toTeX) <- gsub("EC", "Eff. Corr.", rownames(toTeX))
rownames(toTeX) <- gsub("VA", " \\\\& Var. Adj.", rownames(toTeX))

grps <-
  sapply(grps.list, function(x) ifelse(length(x)==3,
                                       paste(x[1], "vs", x[2], "+", x[3]),
                                       paste(x[1], "vs", x[2])))

caption.txt <- "{\\bf CIC data: Comparison of different methods for
        estimating the $\\Delta\\Delta\\textrm{C}_q$-value.}
        In the naive method only undiluted data was used and
        a simple $t$-test performed. LMEM signifies the regular
        $\\Delta\\Delta\\textrm{C}_q$ method using a linaer
        mixed effects model without efficientcy correction.
        Eff.\\ corr.\\ denotes use of the plugin-estimator.
        Var.\\ adj.\\ denotes that the efficientcy correction
        was variance adjusted."
w <- latex(toTeX,
           file    = "../output/Table1.tex",
           title   = "",
           label   = "table:cic",
           caption = caption.txt,
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")
