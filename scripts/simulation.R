
################################################################################
# Simulation example and experiment
# Written by
#   Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted
# Last revision: 10th of June, 2015
################################################################################

#
# Function to perform a simulation
#

SimFunc <- function(nd, ns, n.boots = 151, alpha = 0.05) {

  data <- structure(vector("list", 2), names = c("H0", "HA"))
  for (i in 1:2) {
    data[[i]] <-
      SimqPCRData(std.curve = TRUE, mu.tgt = 25, mu.ref = 30,
                  n.samples = ns, n.replicates = 1, n.dilutions = nd,
                  tech.sd = 1, sample.sd = 1,
                  alpha.tgt = 0.80, alpha.ref = 0.95,
                  ddcq = ifelse(i == 1, 0, 10/9))
  }

  # Aggregate replicates
  w <- FALSE
  qfit0 <- qPCRfit(data$H0, weighted = w)
  qfitA <- qPCRfit(data$HA, weighted = w)

  res <- as.data.frame(
    rbind(
      # Under the null hypothesis
      DDCq.test(data$H0, method = "LMM", eff.cor = FALSE, var.adj = FALSE),
      DDCq(qfit0, eff.cor = TRUE, var.adj = FALSE, alpha = alpha, weighted = w),
      DDCq(qfit0, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w),
      # DDCq(qfit0, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w,
           # var.type = "montecarlo"),
      bs0 <- DDCq.test(data$H0, method =  "Bootstrap", n.boots = n.boots,
                       weighted = w),
      # Under the alternative
      DDCq.test(data$HA, method = "LMM", eff.cor = FALSE, var.adj = FALSE),
      DDCq(qfitA, eff.cor = TRUE, var.adj = FALSE, alpha = alpha, weighted = w),
      DDCq(qfitA, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w),
      # DDCq(qfitA, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w,
           # var.type = "montecarlo"),
      bs1 <- DDCq.test(data$HA, method =  "Bootstrap", n.boots = n.boots,
                       weighted = w)
    ))

  # Sanity checks:
  i <- c(2,6)
  stopifnot(all.equal(res$Estimate[i], res$Estimate[i + 1]))
  if (!all(res$"Std. Error"[i] <= res$"Std. Error"[i + 1])) {
    stop("The standard error in EC+VA is not increased!")
  }
  if (!all(res[i, 5] <= res[i + 1, 5])) {
    stop("The p-values in EC+VA has not increased!")
  }

  rownames(res) <-
    paste0(rep(c("H0:", "H1:"), each = nrow(res)/2),
           rep(c("LMM", "LMM.EC", "LMM.EC.VA1", "LMM.boot"), 2))

  res <- as.matrix(res)
  attr(res, "bootstrapDist") <-
    data.frame(bs0  = attributes(bs0)$extra["Estimate", ],
               bs1  = attributes(bs1)$extra["Estimate", ])
  attr(res, "bs0.warnings") <-  attr(bs0, "warnings")
  attr(res, "bs1.warnings") <-  attr(bs1, "warnings")
  return(res)
}

#
# Perform simulation example
#  - Simulation for 6 samples, 6 dilutions
#

if (!exists("res.ex") || recompute) {

  wrapperSimFunc <- function(seed) {
    set.seed(seed)
    return(as.matrix(SimFunc(nd = 6, ns = 6, n.boots = 201, alpha = 0.05)))
  }

  sfInit(parallel, cpus = n.cpus)
  sfLibrary(lme4)
  sfLibrary(nlme)
  sfExport("wrapperSimFunc", "SimFunc", list = export)

  set.seed(2028674731)  # "Meta-seed": Seed for the seeds

  seeds <- sample.int(2^31, n.sims) # seed for each simulation

  # Do the computation
  res.ex <- sfClusterApplyLB(seeds, wrapperSimFunc)

  sfStop()

  if (any(sapply(res.ex, is.null))) {
    warning("some entires of res.ex are NULL")
  }

  names(res.ex) <- paste0("Sim", seq_along(res.ex))
  resave(res.ex, file = save.file)
}

res.ex <- simplify2array(res.ex)


#
# Perform simulation study
#  - Simulations for combinations of samples and dilutions
#

samples <- c(4, 8)
dilutions <- c(4, 8)

if (!exists("sim.results") || recompute) {

  wrapperSimFunc2 <- function(seed) {
    set.seed(seed)
    return(SimFunc(nd = nd, ns = ns, n.boots = 201, alpha = 0.05))
  }
  st <- proc.time()

  sfInit(parallel, cpus = n.cpus)
  sfLibrary(lme4)
  sfLibrary(nlme)
  sfExport("wrapperSimFunc2", "SimFunc", list = export)

  sim.results        <- vector("list", length(dilutions))
  names(sim.results) <- paste0("n.dilutions", dilutions)
  for (i in seq_along(dilutions)) {
    samp        <- vector("list", length(samples))
    names(samp) <- paste0("n.samples", samples)
    for (j in seq_along(samples)) {

      ns <- samples[j]
      nd <- dilutions[i]
      sfExport("nd", "ns")

      set.seed(840896246)  # "Meta-seed": Seed for the seeds
      seeds <- sample.int(2^31, n.sims) # Seed for each simulation

      res <- sfClusterApplyLB(seeds, wrapperSimFunc2)

      names(res) <- paste0("Sim", seq_along(res))
      samp[[j]] <- simplify2array(res)

      cat(sprintf("dil = %-2d, samp = %-2d, ellapsed: %d mins.\n",
                  i, j, round((proc.time() - st)[3]) %/% 60))
    }
    sim.results[[i]]  <- samp
  }

  sfStop()

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
  hyp <- factor(ifelse(x[,"hypothesis"], "H1", "H0"), levels = c("H1", "H0"))
  sig <- factor(x[,"p"] < p.cut, levels = c(TRUE, FALSE))
  levels(sig) <- sprintf(c("$p < %.2f$", "$p \\geq %.2f$"), p.cut)
  return(table("Truth" = hyp, "Decision" = sig))

}


#
# Create LaTeX table for simulation EXAMPLE
#

p.cut <- 0.05
est <- c("LMM.EC", "LMM.EC.VA1", "LMM.boot")
getr <- paste(rep(c("H0", "H1"), length(est)), rep(est, each = 2), sep = ":")
significant <- res.ex[getr, "Pr(>|t|)", ] < p.cut
ex.tab <- rbind(rowSums(!significant),
                rowSums(significant))
colnames(ex.tab) <- gsub("H0:.+","$H_0$",gsub("H1:.+","$H_A$",colnames(ex.tab)))
rownames(ex.tab) <- sprintf(c("$p \\geq %.2f$", "$p < %.2f$"), p.cut)
ex.cgroup <- c("EC", "EC\\&VA1", "Bootstr.")
tmp.caption <- "Contingency tables for the different estimators for at
  5 \\% $p$-value threshold. The used estimators are the LMM
  with efficiency correction (EC), the LMM
  with EC and variance adjustment (EC\\&VA), and the bootstrapped LMM approach."

w <- latex(ex.tab, file = "../output/Table4.tex", title = "",
           cgroup = ex.cgroup,
           rgroup = "$p$-values",
           caption = tmp.caption,
           label = "tab:simexample")

#
# To use in the knitr document
#

get <- seq(1, ncol(ex.tab) - 1, by = 2)
ex.fpr <- ex.tab["$p < 0.05$", get]/colSums(ex.tab[,get])
ex.tpr <- ex.tab["$p < 0.05$", get + 1]/colSums(ex.tab[,get + 1])
names(ex.fpr) <- names(ex.tpr) <- ex.cgroup

#
# Plot the results
#

setEPS()
postscript("../output/fig4.eps", width = 2*7/1.5, height = 2*7/1.5)

par(mar = c(0, 0, 0, 0), mfrow = c(2, 2), oma = c(5,5.5,2,4), xpd = TRUE)
methods <- c("LMM.EC", "LMM.EC.VA", "LMM.boot")
nm <- length(methods)

h <- 1
for (i in seq_along(dilutions)) {
  for (j in seq_along(samples)) {
    dil <- dilutions[i]
    samp <- samples[j]
    dat <- sim.results[[i]][[j]]

    # Organize data for i and j

    p.cuts <- c(0.01, 0.05, 0.1)
    fpr <- tpr <- as.data.frame(matrix(NA, length(methods)*length(p.cuts), 5))
    names(fpr) <- names(tpr) <- c("p.cut", "est", "rate", "upper", "lower")
    k <- 1
    for (p.cut in p.cuts) {
      for (method in methods) {
        twoByTwo <- t(get.2by2.table(dat, estimator = method, p.cut))

        summ.stats <- summary(epi.tests(twoByTwo))

        fpr[k, ] <- c(p.cut, method, 1 - summ.stats["sp", ])
        tpr[k, ] <- c(p.cut, method, summ.stats["se", ])


        # Check
        dd <- as.data.frame(get.p.info(dat, method), row.names = FALSE)
        dd$sig <- dd$p < p.cut
        dd <- aggregate(sig ~ hypothesis, mean, data = dd)
        stopifnot(all.equal(dd$sig[dd$hypothesis == 0], fpr[k, 3]))
        stopifnot(all.equal(dd$sig[dd$hypothesis == 1], tpr[k, 3]))
        k <- k + 1
      }
    }

    # FPR
    x.fpr <- 1:nrow(fpr) - 0.01
    plot.default(x.fpr, type = "n", axes = FALSE, ylab = "", xlab = "", ylim = 0:1)
    rect(0.5, 0, 1*nm + .5, 1, col = "grey99", border = NA, xpd = TRUE)
    rect(1*nm + .5, 0, 2*nm + .5, 1, col = "grey95", border = NA, xpd = TRUE)
    rect(2*nm + .5, 0, 3*nm + .5, 1, col = "grey90", border = NA, xpd = TRUE)

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

    # Panel letter
    mtext(LETTERS[h], side = 3, adj = 0, cex = 1.1, line = -.5, font = 2)

    # Axes
    if (j == 1) {
      axis(2, las = 2, at = p.cuts, col = "red", col.axis = "red")
      tks <- as.character(axTicks(2))
      tks[tks == "0"] <- ""
      axis(2, las = 2, at = axTicks(2),  label = tks)

      mtext("FPR             ", side = 2, line = 2.5)
      mtext("             TPR", side = 2, line = 2.5, col = col)
      mtext(sprintf("%d dilutions", dil), side = 2, font = 2, line = 4)
    } else {
      axis(4, las = 2, at = p.cuts, col = "red", col.axis = "red")
      tks <- as.character(axTicks(4))
      tks[tks == "0"] <- ""
      axis(4, las = 2, at = axTicks(4),  label = tks)

    }
    if (i == 1) {
      mtext(sprintf("%d samples", samp), side = 3, font = 2, line = 0)
    }

    lab <- gsub("0.[0-9]+ : ", "", gsub("LMM\\.", "", tpr$est))
    lab <- gsub("\\.", "&", lab)
    lab <- gsub("boot", "Boostrap", lab)
    lab <- gsub("VA", "VA1", lab)
    for (d in 1:3) {
      ind <- nm*(d - 1) + 1:nm
      if (i != 1) {
        axis(1, at = ind, labels = lab[ind], las = 2, font = 2)
      }
      # Significance thresholds
      tw <- 0.2
      segments(ind[1] - tw, p.cuts[d], rev(ind)[1] + tw, col = "red", lwd = 1.5)
    }

    # title(main = sprintf("n.samples = %d   n.dilutions = %d", samp, dil))
    h <- h + 1
  }
}
dev.off()


#
# Plot more of the results
#

get.est.info <- function(data, estimator = "LMM.EC") {
  x0 <- cbind(hypothesis = 0, est = data[paste0("H0:", estimator), 1, ])
  x1 <- cbind(hypothesis = 1, est = data[paste0("H1:", estimator), 1, ])
  return(data.frame(x0, x1))
}

setEPS()
postscript("../output/figS1_simCIs.eps", width = 1.5*7, height = 1*7)

par(mfcol = c(6,4), mar = c(0,0,0,0) + 0.2, oma = c(0,4,2,0) + 0.1)
for (i in seq_along(dilutions)) {
  for (j in seq_along(samples)) {
    dil <- dilutions[i]
    samp <- samples[j]
    dat <- sim.results[[i]][[j]]

    # Organize data for i and j
    methods <- c("LMM", "LMM.EC", "LMM.EC.VA")
    p.cuts <- c(0.01, 0.05, 0.1)
    g <- 2

    for (h in seq_along(methods)) {
      for (hyp in c("H0", "H1")) {
        get   <- paste0(hyp, ":", methods[h])
        df    <- dat[get, "df", 1]
        est   <- dat[get, 1, ]
        sd    <- dat[get, 2, ]
        lower <- est + qt(p.cuts[g]/2, df)*sd
        upper <- est - qt(p.cuts[g]/2, df)*sd
        p     <- dat[get, 2, ]
        sig   <- (lower > 0 | upper < 0)
        sig2 <- dat[get, 5, ] < p.cuts[g]
        stopifnot(all(sig == sig2))
        # tp <- (lower < 10/9 & 10/9 < upper) & hyp == "H1"


        plot(est, ylim = c(-4, 5), pch = 16, cex = 0.3, axes = FALSE,
             xlab = "", ylab = "")

        if (hyp == "H0") {
          col <- ifelse(sig, "red", "green")
          legend <- sprintf("FPR = %.3f", mean(sig))
        } else {
          col <- ifelse(sig, "green", "red")
          legend <- sprintf("TPR = %.3f", mean(sig))
        }


        if (j == 1 & i == 1) axis(2)
        segments(seq_along(est), lower, y1 = upper, lwd = 2,
                 col = col)
        abline(h = c(0, 10/9), lty = 1:2)

        if (h == 1 & hyp == "H0") {
          mtext(sprintf("%d dilutions, %d samples", dilutions[i], samples[j]),
                side = 3, outer = FALSE, adj = 0.5, xpd = TRUE, font = 2,
                cex = 0.7)
        }
        if (i == 1 && j == 1) {
          m <- methods[h]
          m <- gsub("\\.", "+", gsub("LMM\\.", "", m))
          mtext(sprintf("%s, %s", m, hyp), line = 3, cex = 0.7,
                side = 2, outer = FALSE, adj = 0.5, xpd = TRUE, font = 2)
        }

        legend("topright", bty = "n", legend = legend)


      }
    }
  }
}
title(paste(100*(1 - p.cuts[g]), "% CIs"), outer = TRUE)
dev.off()





#
# As function of dilutions
#

SimFuncAlt <- function(nd, ns, alpha = 0.05) {

  data <- structure(vector("list", 2), names = c("H0", "HA"))
  for (i in 1:2) {
    data[[i]] <-
      SimqPCRData(std.curve = TRUE, mu.tgt = 25, mu.ref = 30,
                  n.samples = ns, n.replicates = 1, n.dilutions = nd,
                  tech.sd = 1, sample.sd = 1,
                  alpha.tgt = 0.80, alpha.ref = 0.95,
                  ddcq = ifelse(i == 1, 0, 10/9))
  }

  # Aggregate replicates
  w <- FALSE
  qfit0 <- qPCRfit(data$H0, weighted = w)
  qfitA <- qPCRfit(data$HA, weighted = w)

  res <- as.data.frame(
    rbind(
      # Under the null hypothesis
      DDCq.test(data$H0, method = "LMM", eff.cor = FALSE, var.adj = FALSE),
      DDCq(qfit0, eff.cor = TRUE, var.adj = FALSE, alpha = alpha, weighted = w),
      DDCq(qfit0, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w),
      # Under the alternative
      DDCq.test(data$HA, method = "LMM", eff.cor = FALSE, var.adj = FALSE),
      DDCq(qfitA, eff.cor = TRUE, var.adj = FALSE, alpha = alpha, weighted = w),
      DDCq(qfitA, eff.cor = TRUE, var.adj = TRUE, alpha = alpha, weighted = w)
    ))

  # Sanity checks:
  i <- c(2,5)
  stopifnot(all.equal(res$Estimate[i], res$Estimate[i + 1]))
  if (!all(res$"Std. Error"[i] <= res$"Std. Error"[i + 1])) {
    stop("The standard error in EC+VA is not increased!")
  }
  if (!all(res[i, 5] <= res[i + 1, 5])) {
    stop("The p-values in EC+VA has not increased!")
  }

  rownames(res) <-
    paste0(rep(c("H0:", "H1:"), each = nrow(res)/2),
           rep(c("LMM", "LMM.EC", "LMM.EC.VA1"), 2))

  res <- as.matrix(res)
  return(res)
}


dils <- 4:9
samps <- c(4, 7, 10)
if (!exists("sim.results.alt") || recompute) {

  wrapperSimFuncTmp <- function(seed) {
    set.seed(seed)
    return(SimFuncAlt(nd = nd, ns = ns, alpha = 0.05))
  }

  st <- proc.time()

  sfInit(parallel, cpus = n.cpus)
  sfLibrary(lme4)
  sfLibrary(nlme)
  sfExport("wrapperSimFuncTmp", "SimFuncAlt", list = export)

  sim.results.alt        <- vector("list", length(dils))
  names(sim.results.alt) <- paste0("nd", dils)

  for (i in seq_along(dils)) {
    sim.results.alt[[i]] <- vector("list", length(samps))
    names(sim.results.alt[[i]]) <- paste0("ns", samps)

    for (j in seq_along(samps)) {
      ns <- samps[j]
      nd <- dils[i]
      sfExport("nd", "ns")

      set.seed(1437960867)  # "Meta-seed": Seed for the seeds
      seeds <- sample.int(2^31, n.sims) # Seed for each simulation

      sim.results.alt[[i]][[j]] <- sfClusterApplyLB(seeds, wrapperSimFuncTmp)

      cat(sprintf("dil = %-2d, samp = %-2d, ellapsed: %d mins.\n",
                  nd, ns, round((proc.time() - st)[3]) %/% 60))
    }
  }

  sfStop()

  run.time <- proc.time() - st
  attr(sim.results.alt, "time") <- run.time[3]
  resave(sim.results.alt,  file = save.file)
}


# Compute averages across simulations
get.mean <- function(x) {
  return(Reduce("+", x)/length(x))
}

sim.res.alt <- lapply(sim.results.alt, function(x) lapply(x, get.mean))

# Get std errors
get.est.sd <- function(x) {
  return(x[-c(1,4), "Std. Error"])
}

# Format results
sim.res.alt <- lapply(sim.res.alt, function(y) sapply(y, get.est.sd))
sim.res.alt <- simplify2array(sim.res.alt)


# Plot results
setEPS()
postscript("../output/fig5.eps", width = 2*7/1.5, height = 7/1.5)
par(mfrow = c(1,2), mar = c(4, 4, 2, 0) + 0.1)
for (j in 1:2) {

  H <- ifelse(j == 1, "H0", "H1")

  sim.tmp <- sim.res.alt[paste0(H, ":LMM.EC"), , ]
  sim.tmp.va <- sim.res.alt[paste0(H, ":LMM.EC.VA1"), , ]

  samples   <- as.numeric(gsub("ns([0-9]+)",   "\\1", rownames(sim.tmp)))
  dilutions <- as.numeric(gsub("nd([0-9]+)", "\\1", colnames(sim.tmp)))
  cols <- c("steelblue", "darkgrey", "coral")

  ylab <- ifelse(j == 1, expression(paste("Mean SE of ", Delta*Delta*C[q])), "")
  plot(1, type = "n", main = "", ylim = c(0.5, 2), axes = FALSE, log = "y",
       xlim = range(dilutions), xlab = "Dilutions",
       ylab = ylab)
  grid()

  if (j == 1) par(mar = c(4, 3, 2, 0) + 0.1)

  axis(2)
  axis(1, at = dilutions)
  box()

  for (i in seq_len(nrow(sim.tmp))) {
    lines(dilutions, sim.tmp[i, ], col = cols[i],
          lty = 1, lwd = 2, type = "b", pch = 16)
    lines(dilutions, sim.tmp.va[i, ], col = cols[i],
          lty = 2, lwd = 2, type = "b", pch = 1)
  }
  legend("topright", col = cols, lwd = 4,
         legend = gsub("ns([0-9]+)", "\\1 samples", rownames(sim.tmp)),
         lty = 1, bty = "n", inset = 0.025,
         title = parse(text = ifelse(H == "H0", "H[0]", "H[A]")))
  legend("bottomleft", lty = c(1,2), pch = c(16, 1), legend = c("EC", "EC&VA1"),
         bty = "n", inset = 0.025, lwd = 2)

  mtext(LETTERS[j], side = 3, adj = 0, cex = 1.1, line = 1, font = 2)
}
dev.off()




