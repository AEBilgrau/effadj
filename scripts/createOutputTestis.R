
################################################################################
# Analysis of testis data                                                      #
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
  panel.text(2, 40, labels = sprintf("se = %.3f", sqrt(vcov(aa)[2,2])),
             cex = 0.7)
}


# Trellis plot of all data
red.colours <- c("#F82A1C", "#F59691", "#DF6F38", "#AA5455")
stopifnot(nlevels(testis$geneName) == length(red.colours))

testis.data <- subset(testis, sampleType != "Standard")
fig2a <-
  dotplot(sampleName ~ Cq | geneType:sampleType,
          groups = geneName,
          col = red.colours,
          data = testis.data,
          main = "",
          xlab = expression(C[q]),
          pch = seq_along(red.colours),
          key = list(text = list(title = "A"),
                     corner = c(-0.1,1.1),
                     cex = 1.5, font = "bold"))

testis.std <- subset(testis, sampleType == "Standard")
fig2b <-
  xyplot(Cq ~ l2con | geneName,
         groups = geneName,
         col = red.colours,
         data = testis.std,
         panel = mypanel,
         xlab = "Dilution step",#as.expression(bquote(-log[2]*N["0,i,j,k"])),
         ylab = expression(C[q]),
         main = "",
         pch = seq_along(red.colours),
         key = list(text = list(title = "B"),
                    corner = c(-0.1,1.1),
                    cex = 1.5, font = "bold"))

setEPS()
postscript("../output/fig2.eps", width = 1.5*7, height = 0.5*7, fonts = "serif")
  trellis.par.set(strip.background = list(col = "lightgrey"))
  print(fig2a, position = c(0, 0, 0.5, 1), more = TRUE)
  print(fig2b, position = c(0.5, 0, 1, 1))
dev.off()

rm(testis.data, testis.std)

#
# Analysis
#

#
# We wish to test
#   mir127 vs rnu6b,      mir127 vs rnu24,     mir127 vs rnu6b + rnu24 (omitted)
#   mir143 vs rnu6b,      mir143 vs rnu24,     mir143 vs rnu6b + rnu24 (omitted)
#

grps.list <- list(c("mir127", "rnu6b"),
                  c("mir127", "rnu24"),
                  c("mir143", "rnu6b"),
                  c("mir143", "rnu24"))

# Do boostrap
if (!exists("testis.boot") || recompute) {
  message("Testis boostrap")
  testis.boot <- list()
  for (i in 1:length(grps.list)) {
    # Subset data
    testis.tmp <- as.data.qPCR(subset(testis, geneName %in% grps.list[[i]]))

    # Compute bootstrap estimate and results
    testis.boot[[i]] <- bootstrapEstimate(testis.tmp, n.boots = n.boots,
                                          weighted = TRUE, alpha = 0.05)

    message(sprintf("i = %d", i))
  }
  resave(testis.boot, file = save.file)
}


# Combine results
sink("../output/fit_summary_testis.txt")
pdf("../output/modelchecks_testis.pdf", onefile = TRUE)
toTeX <- NULL
we <- TRUE
for (i in 1:length(grps.list)) {
  # Subset data
  testis.tmp <- subset(testis, geneName %in% grps.list[[i]])

  # Print fit information
  m <- paste(grps.list[[i]], collapse = " vs. ")
  cat("\n\n\n\n\n===", m, "===\n")
  print(summary(fit <- qPCRfit(as.data.qPCR(testis.tmp), weighted = we)))
  print(plot(fit, col = testis.tmp$sampleName, pch = testis.tmp$sampleType,
             main = m))

  # Create results for table
  results <- rbind(
    "t-test" = DDCq.test(testis.tmp, method = "N"),
    "LMEM"   = DDCq.test(testis.tmp, method = "LMM",
                         eff.cor = F, weighted = we),
    "EC"     = DDCq.test(testis.tmp, method = "LMM",
                         eff.cor = T, var.adj = F, weighted = we),
    "ECVA1"  = DDCq.test(testis.tmp, method = "LMM",
                         eff.cor = T, var.adj = T, weighted = we),
    "ECVA2"  = DDCq.test(testis.tmp, method = "LMM",
                         eff.cor = T, var.adj = T,
                         var.type = "montecarlo", weighted = we),
    "Bootstrap" = testis.boot[[i]]
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
  sapply(grps.list, function(x) ifelse(length(x) == 3,
                                       paste(x[1], "vs", x[2], "+", x[3]),
                                       paste(x[1], "vs", x[2])))

caption.txt <- "Testis data: Method comparison for estimating the
  $\\ddcq$-value.
  EC denotes use of the plugin-estimator.
  VA denotes that the efficiency correction was variance adjusted using the
  delta method (1) or Monte Carlo integration (2).
  Bootstrap shows the mean and standard deviation of %d
  bootstrap samples using EC estimate. The last two columns show the $95%s$
  lower and upper confidence interval limits."

toTeX <- toTeX[!grepl("LMM|t-test", rownames(toTeX)), ]

w <- latex(toTeX,
           file    = "../output/Table2.tex",
           title   = "",
           label   = "table:tesits",
           caption = sprintf(caption.txt, length(testis.boot), "\\%"),
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")





