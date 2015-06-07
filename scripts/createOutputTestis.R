
################################################################################
# Analyse testis and create outpt                                              #
# Written by Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted      #
# Last revision: 7th of June, 2015                                             #
################################################################################

#
# Initalization
#

# Removing mir143 and rnu24
#testis <- testis[!(testis$geneName %in% c("mir143", "rnu24")), ]

# Trellis plot of all data miR-127, miR-143, miR-6b, miR-24 data
pdf("../output/Figure3.testis.pdf", width = 7*2, height = 4)
  trellis.par.set(strip.background=list(col="lightgrey"))

  fig3a <- dotplot(sampleName ~ Cq | geneType:sampleType,
                   data = testis[testis$sampleType != "Standard", ],
                   main = "", col = "tomato",
                   key = list(text=list(title="A"),
                              corner = c(-0.1,1.1),
                              cex = 1.5, font = "bold"))

  fig3b <- dotplot(l2con ~ Cq | geneType:geneName,
                   scales = list(y = list(lab = sort(unique(testis$l2con)))),
                   data = testis[testis$sampleType == "Standard",],
                   ylab = as.expression(bquote(-log[2]*N["0,i,j,k"])),
                   main = "", col = "tomato",
                   key = list(text=list(title="B"),
                              corner = c(-0.1,1.1),
                              cex = 1.5, font = "bold"))

  print(fig3a, position=c(0, 0, 0.5, 1), more=TRUE)
  print(fig3b, position=c(0.5, 0, 1, 1))
dev.off()


#
# Bootstrapped confidence bounds and p-values
#

new.testis <- aggregate(Cq ~ sampleName + geneName + sampleType +
                          geneType + l2con + copyNumber,
                        data = testis, FUN = mean)

new.testis$replicate          <- 1
class(new.testis)             <- c("data.qPCR", "data.frame")
attr(new.testis, "std.curve") <- TRUE

start.sample   <- 3
n.samples      <- 5
start.dilution <- 3
n.dilutions    <- 5
end.sample <- start.sample + n.samples - 1

if (!exists("testis.boot") || recompute) {
  testis.boot  <- Bootstrap.qPCR(new.testis,
                                 start.sample   = start.sample,
                                 n.samples      = n.samples,
                                 start.dilution = start.dilution,
                                 n.dilutions    = n.dilutions,
                                 n.resamp       = 1000)
  resave(testis.boot, file = save.file)
}

# Plotting results
pdf("../output/testis.pcurve.pdf")

  plot(1, type = "n",
       xlab = "Number of samples per group",
       xlim = c(start.sample, end.sample),
       ylab = "log P-value",
       ylim = c(0.01, 1), log = "")
  grid()
  lines(c(start.sample, end.sample), c(0.05, 0.05), lty = 2)

  for (i in 1:4) {
    lines(start.sample:end.sample, testis.boot[[i]][4,],
          type = "b", col = i, lwd = 1.5)
    lines(start.sample:end.sample, testis.boot[[i]][5,],
          type = "b", col = i, lwd = 1.5, lty = 2)
    lines(start.sample:end.sample, testis.boot[[i]][6,],
          type = "b", col = i, lwd = 1.5, lty = 2)
  }

  legend("topright",
         legend = c("4 step dilution", "5 step dilution",
                    "6 step dilution", "7 step dilution"),
         lty    = rep(1, 4),
         lwd    = rep(1.5, 4),
         col    = 1:4,
         bty    = "n",
         title  = "")

dev.off()


pdf("../output/testis.estcurve.pdf")

  plot(1, type = "n",
       xlab = "Number of samples per group",
       xlim = c(start.sample, start.sample + n.samples - 1),
       ylab = "DDCq",
       ylim = c(0,4))
  grid()

  for(i in 1:4){
    lines(start.sample:end.sample,testis.boot[[i]][1,],
          type = "b", col = i, lwd = 1.5)
    lines(start.sample:end.sample,testis.boot[[i]][2,],
          type = "b", col = i, lwd = 1.5, lty = 2)
    lines(start.sample:end.sample,testis.boot[[i]][3,],
          type = "b", col = i, lwd = 1.5, lty = 2)
  }

  legend("bottomright",
         legend = c("4 step dilution", "5 step dilution",
                    "6 step dilution", "7 step dilution"),
         lty    = rep(1, 4),
         lwd    = rep(1.5, 4),
         col    = 1:4,
         bty    = "n",
         title  = "")

dev.off()

#
# Analysis
#


# We wish to test
# mir127 vs rnu6b
# mir127 vs rnu24
# mir127 vs rnu6b + rnu24

grps.list <- list(c("mir127", "rnu6b"),
                  c("mir127", "rnu24"),
                  c("mir127", "rnu6b", "rnu24"),
                  c("mir143", "rnu6b"),
                  c("mir143", "rnu24"),
                  c("mir143", "rnu6b", "rnu24"))

toTeX <- NULL
for (i in 1:length(grps.list)) {

  genes <- grps.list[[i]]
  testis.tmp <- as.data.qPCR(testis[testis$geneName %in% genes, ])

  results <-
    rbind("t-test" = DDCq.test(testis.tmp, method = "N"),
          "LMEM"   = DDCq.test(testis.tmp, method = "LMM", eff.cor = FALSE),
          "EC"     = DDCq.test(testis.tmp, method = "LMM",
                               eff.cor = TRUE, var.adj = FALSE),
          "ECVA"   = DDCq.test(testis.tmp, method = "LMM",
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

caption.txt <- "{\\bf Testis data: Comparison of different methods for
estimating the $\\Delta\\Delta\\textrm{C}_q$-value.}
In the naive method only undiluted data was used and
a simple $t$-test performed. LMEM signifies the regular
$\\Delta\\Delta\\textrm{C}_q$ method using a linaer
mixed effects model without efficientcy correction.
Eff.\\ corr.\\ denotes use of the plugin-estimator.
Var.\\ adj.\\ denotes that the efficientcy correction
was variance adjusted."

w <- latex(toTeX,
           file    = "../output/Table2.tex",
           title   = "",
           label   = "table:tesits",
           caption = caption.txt,
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")
