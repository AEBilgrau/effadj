
################################################################################
# Simulation example and experiment
# Written by
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted
# Last revision: 10th of June, 2013
################################################################################

#
# Perform a simulation
#

SimTemp <- function (nd, ns, ddcq = 10/9) {
  # Initalization / parameters of simulation
  mu.tgt         <- 30     # Mean of target gene
  mu.ref         <- 25     # Mean of reference gene
  alpha.tgt      <- 0.75   # Amp. efficiency of target gene
  alpha.ref      <- 0.90   # Amp. efficiency of refence gene
  ddcq           <- 10/9   # True effect size Delta-Delta-C_q
  tech.sd        <- 1      # Technical standard deviation
  sample.sd      <- 1      # Sample standard deviation
  n.replicates   <- 1      # Number of technical replicates
  n.boot.sim     <- 150    # Number of bootstrap samples

  data <- vector("list", 2)
  for (i in 1:2) {
    data[[i]] <-
      SimqPCRData(std.curve = TRUE, mu.tgt = mu.tgt, mu.ref = mu.ref,
                  n.samples = ns, n.replicates = n.replicates, n.dilutions = nd,
                  tech.sd = tech.sd, alpha.tgt = alpha.tgt,
                  alpha.ref = alpha.ref, sample.sd = sample.sd,
                  ddcq = ifelse(i==1, 0, ddcq))
  }
  res <- as.data.frame(
    rbind(DDCq.test(data[[1]], method = "LMM", eff.cor=FALSE, var.adj=FALSE),
          DDCq.test(data[[1]], method = "LMM", eff.cor=TRUE,  var.adj=FALSE),
          DDCq.test(data[[1]], method = "LMM", eff.cor=TRUE,  var.adj=TRUE),
          DDCq.test(data[[1]], method = "Bootstrap", n.boot.sim),
          DDCq.test(data[[2]], method = "LMM", eff.cor=FALSE, var.adj=FALSE),
          DDCq.test(data[[2]], method = "LMM", eff.cor=TRUE,  var.adj=FALSE),
          DDCq.test(data[[2]], method = "LMM", eff.cor=TRUE,  var.adj=TRUE),
          DDCq.test(data[[2]], method = "Bootstrap", n.boot.sim)))

  # Sanity checks:
  stopifnot(all.equal(res$Estimate[c(2,6)], res$Estimate[c(2,6)+1]))
  if (!all(res$"Std. Error"[c(2,6)] <= res$"Std. Error"[c(2,6)+1])) {
    stop("The standard error in EC+VA is not increased!")
  }
  if (!all(res[c(2,6), 5] <= res[c(2,6)+1, 5])) {
    stop("The p-values in EC+VA has not increased!")
  }

  rownames(res) <-
    paste0(rep(c("H0:", "H1:"), each = 4),
           rep(c("LMM", "LMM.EC", "LMM.EC.VA", "LMM.boot"), 2))
  return(res)
}

#
# Perform simulation example
#  - Simulation for 6 samples, 6 dilutions
#

if (!exists("res.ex") || recompute) {
  ex <- SimTemp(nd = 3, ns = 3)  # Just used to get dimnames hereof
  res.ex <- array(NA, c(8, 5, n.sims))
  dimnames(res.ex) <- c(dimnames(ex), list(paste0("Sim", seq_len(n.sims))))
  for (k in seq_len(n.sims)) {
    res.ex[, , k] <- as.matrix(SimTemp(nd = 6, ns = 6))
    cat(sprintf("sim = %d\n", k))
  }
  resave(res.ex, file = save.file)
}

#
# Perform simulation study
#  - Simulations for combinations of samples and dilutions
#

samples <- c(4, 8)
dilutions <- c(4, 8)

if (!exists("sim.results") | recompute) {
  ex <- SimTemp(nd = 3, ns = 3)  # Just used to get dimnames hereof
  set.seed(36)
  st <- proc.time()

  sim.results        <- vector("list", length(dilutions))
  names(sim.results) <- paste0("n.dilutions", dilutions)
  for (i in seq_along(dilutions)) {
    samp        <- vector("list", length(samples))
    names(samp) <- paste0("n.samples", samples)
    for (j in seq_along(samples)) {
      res <- array(NA, c(8, 5, n.sims))
      dimnames(res) <- c(dimnames(ex), list(paste0("Sim", seq_len(n.sims))))

      for (k in seq_len(n.sims)) {
        res[, , k] <- as.matrix(SimTemp(nd = dilutions[i], ns = samples[j]))

        tm <- (proc.time() - st)/60
        cat(sprintf("dil = %-2d, samp = %-2d, sim = %-3d, ellapsed: %d mins.\n",
                    i, j, k, round(tm[3])))
      }
      samp[[j]] <- res
    }
    sim.results[[i]]  <- samp
  }

  run.time <- proc.time() - st
  attr(sim.results, "time") <- run.time[3]
  resave(sim.results,  file = save.file)
}

#
# Plot and write results
#

get.p.info   <- function(subdata, estimator = "LMM.EC") {
  x0 <- cbind(hypothesis = 0, p = subdata[paste0("H0:", estimator), 5, ])
  x1 <- cbind(hypothesis = 1, p = subdata[paste0("H1:", estimator), 5, ])
  return(rbind(x0, x1))
}

# Sanity check (again)
for (i in seq_along(dilutions)) {
  for (j in seq_along(samples)) {
    A <- get.p.info(sim.results[[i]][[j]], estimator = "LMM.EC")
    B <- get.p.info(sim.results[[i]][[j]], estimator = "LMM.EC.VA")
    stopifnot(all(A[,2] < B[,2]))
  }
}

get.2by2.table <-  function(subdata, estimator = "LMM.EC", p.cut = 0.05) {
  x <- get.p.info(subdata, estimator)
  hyp <- ifelse(x[,"hypothesis"], "H1", "H0")
  sig <- x[,"p"] < p.cut
  dec <- sprintf(ifelse(sig, "$p < %.2f$", "$p \\geq %.2f$"), p.cut)
  return(table("Truth" = hyp, "Decision" = dec))
}

get.performance <- function(subdata, estimator){
  x <- get.p.info(subdata, estimator)
  pred <- prediction(predictions = x[, "p"], labels = x[, "hypothesis"])
  tpr <- performance(pred, "tpr", "cutoff")
  fpr <- performance(pred, "fpr", "cutoff")
  return(list("tpr" = tpr, "fpr" = fpr))
}



# Create LaTeX table for simulation EXAMPLE
tab <- matrix(NA, 2, 0)
est <- c("LMM", "LMM.EC", "LMM.EC.VA", "LMM.boot")
for (nm in est) {
  subtab <- t(get.2by2.table(res.ex, estimator = nm, p.cut = 0.05))[2:1, ]
  tab <- cbind(tab, subtab)
}
colnames(tab) <- gsub("H0", "$H_0$", gsub("H1", "$H_A$", colnames(tab)))
tmp.caption <- "Contingency tables for the different estimators for a $p$-value
  cutoff of 0.05."
w <- latex(tab, file = "../output/Table3.tex", title = "",
           cgroup = c("LMM", "EC", "EC.VA", "Bootstrap"),
           rgroup = "Decision",
           caption = tmp.caption,
           label = "tab:simexample")


#
# Plot the results
#

setEPS()
postscript("../output/fig3.eps", width = 2*7/1.5, height = 2*7/1.5)

par(mar = c(0, 0, 0, 0) + 0.5, mfrow = c(2, 2), oma = c(4.5,5.5,2,0), xpd=TRUE)

for (i in seq_along(dilutions)) {
  for (j in seq_along(samples)) {
    dil <- dilutions[i]
    samp <- samples[j]
    dat <- sim.results[[1]][[2]]

    # Organize data for i and j
    methods <- c("LMM", "LMM.EC", "LMM.EC.VA", "LMM.boot")
    p.cuts <- c(0.01, 0.05, 0.1)
    fpr <- tpr <- as.data.frame(matrix(NA, length(methods)*length(p.cuts), 5))
    names(fpr) <- names(tpr) <- c("p.cut", "est", "rate", "upper", "lower")
    k <- 1
    for (p.cut in p.cuts) {
      for (method in methods) {
        twoByTwo <- get.2by2.table(dat, estimator = method, p.cut)
        summ.stats <- summary(epi.tests(t(twoByTwo)[2:1, 2:1]))
        fpr[k, ] <- c(p.cut, method, 1 - summ.stats["sp", ])
        tpr[k, ] <- c(p.cut, method, summ.stats["se", ])
        k <- k + 1
      }
    }

    # FPR
    x.fpr <- 1:nrow(fpr) - 0.01
    plot.default(x.fpr, type = "n", axes = FALSE, ylab="", xlab="", ylim = 0:1)
    segments(x.fpr, fpr$upper, y1 = fpr$lower)
    points(x.fpr, fpr$upper, pch = "-")
    points(x.fpr, fpr$lower, pch = "-")
    points(x.fpr, fpr$rate,  pch = 15:18)

    # TPR
    x.tpr <- 1:nrow(tpr) + 0.02
    col <- "darkgrey"
    segments(x.tpr, tpr$upper, y1 = tpr$lower, lty = 2, col = col)
    points(x.tpr, tpr$upper, pch = "-",   col = col)
    points(x.tpr, tpr$lower, pch = "-",   col = col)
    points(x.tpr, tpr$rate,  pch = 15:18, col = col)

    if (j == 1) {
      axis(2)
      mtext("FPR             ", side = 2, line = 2)
      mtext("             TPR", side = 2, line = 2, col = col)
      mtext(sprintf("%d dilutions", dil), side = 2, font = 2, line = 4)
    }
    if (i == 1) {
      mtext(sprintf("%d samples", samp), side = 3, font = 2, line = 0)
    }

    lab <- gsub("0.[0-9]+ : ", "", gsub("LMM\\.", "", tpr$est))
    for (d in 1:3) {
      ind <- 4*(d-1) + 1:4
      if (i == 1) {lab <- ""}
      axis(1, at = ind, labels = lab[ind], las = 2)
      tw <- 0.1

      # Significance thresholds
      segments(ind[1]-tw, c(0.01, 0.05, 0.1)[d], rev(ind)[1]+tw,
               col = "red", lwd = 2)
    }

    # title(main = sprintf("n.samples = %d   n.dilutions = %d", samp, dil))
  }
}
dev.off()

