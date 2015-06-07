
################################################################################
# Simulations
# By: Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted
# Last revision: 10th of June, 2013
################################################################################

#
# Initalization
#

mu.tgt         <- 30     # Mean of target gene
mu.ref         <- 25     # Mean of reference gene
alpha.tgt      <- 0.75   # Amp. efficiency of target gene
alpha.ref      <- 0.80   # Amp. efficiency of refence gene
ddcq           <- 10/9   # True effect size Delta-Delta-C_q
tech.sd        <- 1      # Technical standard deviation
sample.sd      <- 1.3    # Sample standrad deviation
n.sims         <- 100 #10000   # Number of simulations
n.replicates   <- 1      # Number of technical replicates

start.sample   <- 3      # Smallest number of samples
n.samples      <- 10     # Total number of samples
end.sample     <- start.sample + n.samples - 1 # Max number of samples
start.dilution <- 3      # Smallest number of dilutions
n.dilutions    <- 12     # Total number of dilutions


#
# Power calculations
#

# Power calculation with standard curves
if (!exists("dilution.power.results") | recompute) {
  dilution.power.results <-
    PowerSim(n.sims         = n.sims,
             start.sample   = start.sample,
             n.samples      = n.samples,
             start.dilution = start.dilution,
             n.dilutions    = n.dilutions,
             ddcq           = ddcq,
             alpha.tgt = alpha.tgt, alpha.ref = alpha.ref,
             mu.tgt = mu.tgt, mu.ref = mu.ref, tech.sd = tech.sd,
             n.replicates = n.replicates)
  resave(dilution.power.results, file = save.file)
}

Rprof()
# Power calculations without standard curves
if (!exists("no.dilution.power.results") | recompute) {
  no.dilution.power.results <-
    PowerSim(std.curve      = FALSE,
             n.sims         = n.sims,
             start.sample   = start.sample,
             n.samples      = n.samples,
             start.dilution = 1,        # These are equivalent
             n.dilutions    = 1,        # to std.curve = FALSE
             ddcq           = ddcq,
             alpha.tgt = alpha.tgt, alpha.ref = alpha.ref,
             mu.tgt = mu.tgt, mu.ref = mu.ref, tech.sd = tech.sd,
             n.replicates = n.replicates)
  save(no.dilution.power.results, file = save.file)
}
Rprof(NULL)
summaryRprof()


# Theoretical power curve, with perfect efficiency
t.pow.res <- rep(NA, n.samples)
for(i in 1:length(t.pow.res)){
  t.pow.res[i] <- power.t.test(n           = start.sample + i - 1,
                               delta       = ddcq,
                               sd          = sample.sd,
                               sig.level   = 0.05,
                               power       = NULL,
                               type        = c("two.sample"),
                               alternative = c("two.sided"),
                               strict      = FALSE)$power
}

# Plotting power curves

pow.res     <- dilution.power.results
without.res <- no.dilution.power.results

pdf("output/Figure1.pdf")
  plot(1, type = "n",
       xlab = "Number of samples per group",
       xlim = c(start.sample, end.sample),
       ylab = "Power",
       ylim = c(min(pow.res,t.pow.res, without.res),
                max(pow.res,t.pow.res, without.res)))
  grid()
  for (i in 1:nrow(pow.res)) {
    lines(cenvelope(cbind(start.sample:end.sample, pow.res[i, ])),
          type = "b", col = jet.colors(5)[i], lwd = 1.5, cex = 0.6)
  }
  lines(start.sample:end.sample, without.res[1,],
        type = "b", col = "grey", lwd = 1.5, lty = 2, cex = 0.6)
  lines(start.sample:end.sample, t.pow.res,
        type = "b", col = "black", lwd = 1.5, lty = 2, cex = 0.6)


  legend("bottomright",
         legend = c(rownames(pow.res), "No dilution", "t-test approach"),
         lty = c(rep(1, nrow(pow.res)), 2, 2), lwd = 2, bty = "n",
         col = c(jet.colors(nrow(pow.res)), "grey", "black"), inset = 0.0)
  legend("right", legend = bquote(Delta*Delta*C[q] == .(round(ddcq,3))),
         bty = "n", inset = 0.05)
dev.off()


#
# More simulations
#


SimTemp <- function (nd, ns, ddcq = 10/9) {
  data <- list()
  for (i in 1:2) {
    data[[i]] <-
      SimqPCRData(std.curve = TRUE, mu.tgt = mu.tgt, mu.ref = mu.ref,
                  n.samples = ns, n.replicates = n.replicates, n.dilutions = nd,
                  tech.sd = tech.sd, alpha.tgt = alpha.tgt, alpha.ref = alpha.ref,
                  sample.sd = sample.sd, tech.sd = tech.sd,
                  ddcq = ifelse(i==1, 0, ddcq))
  }
  res <- as.data.frame(
    rbind(DDCq.test(data[[1]], method = "LMM", eff.co =FALSE, var.adj=FALSE),
          DDCq.test(data[[1]], method = "LMM", eff.cor=TRUE,  var.adj=FALSE),
          DDCq.test(data[[1]], method = "LMM", eff.cor=TRUE,  var.adj=TRUE),
          DDCq.test(data[[2]], method = "LMM", eff.cor=FALSE, var.adj=FALSE),
          DDCq.test(data[[2]], method = "LMM", eff.cor=TRUE,  var.adj=FALSE),
          DDCq.test(data[[2]], method = "LMM", eff.cor=TRUE,  var.adj=TRUE)))
  rownames(res) <- paste(rep(c("H0:", "H1:"), each = 3),
                         rep(c("LMM", "LMM.EC", "LMM.EC.VA"), 2), sep = "")
  return(res)
}

# Parameters controling the number of samples and dilutions to be simulated under
from.samp <- 4
to.samp   <- 11
from.dil  <- 4
to.dil    <- 12

if (!exists("sim.results") | recompute) {

  st <- proc.time()
  sim.results        <- vector("list", to.dil - from.dil + 1)
  names(sim.results) <- paste("dilution", from.dil:to.dil, sep = "")

  for (i in seq_along(sim.results)) {
    samp        <- vector("list", to.samp - from.samp + 1)
    names(samp) <- paste("sample", from.samp:to.samp, sep = "")

    for (j in seq_along(samp)) {
      res <- array(NA, c(6, 5, n.sims))
      cat("Iteration:")
      for (k in seq_len(n.sims)) {
        res[, , k] <-
          as.matrix(SimTemp((from.dil:to.dil)[i],
                            (from.samp:to.samp)[j]))

        if (k%%100==0) cat(k, " "); flush.console()
        if (k%%1000==0) cat("\n")
      }
      samp[[j]] <- res;
      cat("Samples:", j+from.samp-1, "done.\n")
      cat("Dilutions:", i+from.dil-1, "done.\n"); flush.console()
    }
    sim.results[[i]]  <- samp
  }
  resave(sim.results,  file = save.file)
  run.time <- proc.time() - st
}

#
# Plot results
#

get.EC   <- function(data) {
  cbind(rep(0:1, each = dim(data)[3]), c(data[2,5, ], data[5,5, ]))
}

get.ECVA <- function(data) {
  cbind(rep(0:1, each = dim(data)[3]), c(data[3,5, ], data[6,5, ]))
}


#
# Plotting ROC curves
#

get.roc <- function(x){
  pred <- prediction(-x[,2], x[,1])
  perf <- performance(pred, "tpr", "fpr")
}

get.auc <- function(x) {
  pred <- prediction(-x[,2], x[,1])
  unlist(slot(performance(pred, "auc"), "y.values"))
}

roc.ec <- lapply(sim.results, function(x) {lapply(x, get.EC)})
aucs.ec <- sapply(roc.ec, function(x) {sapply(x, get.auc)})
roc.ec <- lapply(roc.ec, function(x) {lapply(x, get.roc)})

roc.ecva <- lapply(sim.results, function(x) {lapply(x, get.ECVA)})
aucs.ecva <- sapply(roc.ecva, function(x) {sapply(x, get.auc)})
roc.ecva <- lapply(roc.ecva, function(x) {lapply(x, get.roc)})


#
# All ROC curves
#

png("output/simulation.roc.curves.png", height = 1.5*7, width = 1.5*7,
    units = "in", res = 300)
par(mfrow=c(3,3))
for (i in 1:length(roc.ec)) {
  for (j in 1:length(roc.ec[[1]])) {

    plot(roc.ec[[i]][[j]], add = ifelse(j==1, FALSE, TRUE),
         col = "red", lwd = 1, lty = 1, main = names(roc.ec)[i])
    plot(roc.ecva[[i]][[j]], add = TRUE,
         col = "blue", lwd = 1, lty = 2)
    if (j==1) {
      abline(0, 1, col = "grey", lty = 2)
      legend("bottomright", inset = 0.02, bty = "n", lty = c(1,2), lwd = 2,
             legend = c("EC", "EC + VA"),
             col = c("red", "blue"))
    }
  }
}
dev.off()

#
# AUC ratios between EC+VA and EC
#

pdf("output/simulation.auc.ratio.pdf")
  heatmap.2(aucs.ecva/aucs.ec, dendrogram  = "non", keysize = 1,
            density.info = "none", main = "AUC(EC+VA) / AUC(EC)",
            trace = "none", Rowv = FALSE, Colv = FALSE)
dev.off()

#
# Threshold against FPR and TPR
#

png("output/simulation.threshold.vs.FPR.TPR.png", width = 14, height = 7,
    units = "in", res = 300)
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  split.screen(c(1,2))
  for (i in 1:length(roc.ec)) {
    for (j in 1:length(roc.ec[[1]])) {

      roc.ec.tmp   <- roc.ec[[i]][[j]]
      roc.ecva.tmp <- roc.ecva[[i]][[j]]

      l <- i == 1 & j == 1
      screen(1)
      plot(1, type = "n", xlim = 0:1, ylim = 0:1,
           xlab = ifelse(l, expression(alpha), ""), ylab = ifelse(l, "FPR", ""),
           axes = l)
      lines(-roc.ec.tmp@alpha.values[[1]], roc.ec.tmp@x.values[[1]],
            type = "s", col = "red", lty = 1)
      lines(-roc.ecva.tmp@alpha.values[[1]], roc.ecva.tmp@x.values[[1]],
            type = "s", col = "blue", lty = 2)

      screen(2)
      plot(1, type = "n", xlim = 0:1, ylim = 0:1,
           xlab = ifelse(l, expression(alpha), ""), ylab = ifelse(l, "TPR", ""),
           axes = l)
      lines(-roc.ec.tmp@alpha.values[[1]], roc.ec.tmp@y.values[[1]],
            type = "s", col = "red", lty = 1)
      lines(-roc.ecva.tmp@alpha.values[[1]], roc.ecva.tmp@y.values[[1]],
            type = "s", col = "blue", lty = 2)

    }
  }
  legend("bottomright", inset = 0.02, bty = "n", lty = c(1,2), lwd = 2,
         legend = c("EC", "EC + VA"),
         col = c("Red", "Blue"))
dev.off()

#
# Plotting more
#

matrixMeans <- function(data) {
  apply(data, c(1,2), mean)
}

sim.tmp <- lapply(sim.results, function(x) lapply(x, matrixMeans))
sim.tmp <- as.data.frame(lapply(sim.tmp, function(x) sapply(x,"[",5,2)))

sim.tmp.va <- lapply(sim.results, function(x) lapply(x, matrixMeans))
sim.tmp.va <- as.data.frame(lapply(sim.tmp.va, function(x) sapply(x,"[",6,2)))

samples   <- as.numeric(gsub("sample([0-9]+)",   "\\1", rownames(sim.tmp)))
dilutions <- as.numeric(gsub("dilution([0-9]+)", "\\1", colnames(sim.tmp)))

pdf("output/simulation.mean.sd.of.ddcq.pdf", height = 7, width = 7)

  plot(1, type = "n", main = "", ylim = c(0.3, 1), axes = FALSE,
       xlim = range(dilutions), xlab = "Dilutions",
       ylab = expression(paste("Mean SD of ", Delta*Delta*C[q])))
  axis(2); axis(1, at = dilutions); grid(); box()
  for (i in 1:nrow(sim.tmp)) {
    lines(dilutions, sim.tmp[i, ],
          col = jet.colours(nrow(sim.tmp))[i])
    lines(dilutions, sim.tmp.va[i, ],
          col = jet.colours(nrow(sim.tmp))[i], lty = 2)
  }
  legend("topright", col = jet.colours(nrow(sim.tmp)),
         legend = gsub("sample([0-9]+)", "\\1 samples", rownames(sim.tmp)),
         lty = 1, bty = "n", inset = 0.025)
  legend("bottomleft", lty = c(1,2), legend = c("EC", "EC & VA"),
         bty = "n", inset = 0.025)
dev.off()


