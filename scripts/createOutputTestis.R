
################################################################################
# Analysis of testis data                                                      #
# Written by:                                                                  #
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted               #
# Last revision: 11th of June, 2015                                            #
################################################################################

#
# Plot raw data
#

# Trellis plot of all data
testis.data <- subset(testis, sampleType != "Standard")
fig2a <- dotplot(sampleName ~ Cq | geneType:sampleType,
                 data = testis.data,
                 main = "", col = "tomato",
                 xlab = expression(C[q]),
                 pch = 16,
                 key = list(text=list(title="A"),
                            corner = c(-0.1,1.1),
                            cex = 1.5, font = "bold"))

testis.std <- subset(testis, sampleType == "Standard")
fig2b <- xyplot(Cq ~ l2con | geneName,
                data = testis.std,
                xlab = as.expression(bquote(-log[2]*N["0,i,j,k"])),
                ylab = expression(C[q]),
                main = "", col = "tomato",
                pch = 16,
                key = list(text=list(title="B"),
                           corner = c(-0.1,1.1),
                           cex = 1.5, font = "bold"))

setEPS()
postscript("../output/fig2.eps", width = 1.5*7, height = 0.5*7)
  trellis.par.set(strip.background=list(col="lightgrey"))
  print(fig2a, position=c(0, 0, 0.5, 1), more = TRUE)
  print(fig2b, position=c(0.5, 0, 1, 1))
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
if (!exists("testis.boot") || !exists("testis.pboot") || recompute) {
  message("Testis boostrap")
  testis.boot <- testis.pboot <- list()
  for (i in 1:length(grps.list)) {
    # Subset data
    testis.tmp <- as.data.qPCR(subset(testis, geneName %in% grps.list[[i]]))

    # Compute bootstrap estimate and results
    testis.boot[[i]] <- bootstrapEstimate(testis.tmp, n.boots = n.boots)
    testis.pboot[[i]]<- parametricBootstrapEstimate(testis.tmp, n.boots=n.boots)

    message(sprintf("i = %d", i))
  }
  resave(testis.boot, testis.pboot, file = save.file)
}


# Combine results
sink("../output/fit_summary_testis.txt")
toTeX <- NULL
for (i in 1:length(grps.list)) {
  # Subset data
  testis.tmp <- subset(testis, geneName %in% grps.list[[i]])

  # Print fit information
  fit <- qPCRfit(as.data.qPCR(cic.tmp))
  cat("\n\n\n\n\n===", paste(grps.list[[i]], collapse = " vs. "),"===\n")
  print(summary(fit))

  # Create results for table
  results <- rbind(
    "t-test" = DDCq.test(testis.tmp, method = "N"),
    "LMEM"   = DDCq.test(testis.tmp, method = "LMM", eff.cor = F),
    "EC"     = DDCq.test(testis.tmp, method = "LMM", eff.cor = T, var.adj = F),
    "ECVA1"  = DDCq.test(testis.tmp, method = "LMM", eff.cor = T, var.adj = T),
    "ECVA2"  = DDCq.test(testis.tmp, method = "LMM", eff.cor = T, var.adj = T,
                         var.type = "montecarlo"),
    "Bootstrap" = testis.boot[[i]],
    "Par.~bootstr." = testis.pboot[[i]]
    )

  toTeX <- rbind(toTeX, results)
}
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
  sapply(grps.list, function(x) ifelse(length(x)==3,
                                       paste(x[1], "vs", x[2], "+", x[3]),
                                       paste(x[1], "vs", x[2])))

caption.txt <- "Testis data: Method comparison for estimating the
  $\\ddcq$-value. $t$-test shows results from a simple unpaired
  $t$-test ignoring dilution data. LMM signifies the regular
  $\\ddcq$ method using a linear mixed effects model ignoring
  dilution.
  EC denotes use of the plugin-estimator.
  VA denotes that the efficiency correction was variance adjusted using the
  delta method (1) or monte carlo integration (2).
  (Par.) bootstrap shows the mean and standard deviation of %d (parametric)
  bootstrap samples using EC estimate. The last two columns shows the $95%s$
  lower and upper confidence interval limits."
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





