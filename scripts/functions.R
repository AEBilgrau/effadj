
################################################################################
# Core and Auxillary functions
# By: Anders Ellern Bilgrau, Steffen Falgreen, and Martin Boegsted
# Last revision: 16th of Feb, 2015
################################################################################

# This script defines the S3 class "data.qPCR". The data.qPCR-class is a
# subclass of data.frames. Therefore, data.qPCR have class
# c("data.qPCR", "data.frame") and the additional logical attribute
# "std.curve" indicating whether standard curve data is present
# in the dataset. In addition, the the object must contain the columns:
#
#   "sampleName"  (factor; wit sample names)
#   "sampleType"  (factor with 2 or 3 levels: "case", "control", "Standard"
#   "geneType"    (factor with 2 levels: "tgt", "ref"
#   "replicate",  (factor; the sample replicate number)
#   "Cq",         (numeric; the Cq (aka Ct) values of the qPCR experiment)
#   "copyNumber", (numeric; the dilution, e.g. 1.0, 0.5, 0.25, etc.)
#   "l2con"       (numeric; minus log2 of the concentration)
#
# Note, a prototype of a data.qPCR instance can be created by running
# SimqPCRData() with (std.curve = TRUE) or without (std.curve = FALSE)
# standard curves.


#
# Simulation of a random data.qPCR dataset
#

SimqPCRData <-
  function(std.curve    = TRUE,      # Logical; simulate standard curves or not
           mu.tgt       = 30,        # Mean of target gene
           mu.ref       = 25,        # Mean of referece gene
           n.samples    = 4,         # Number of samples
           n.replicates = 3,         # Number of replicates
           n.dilutions  = 5,         # Number of dilutions
           tech.sd      = 1,         # Variance of technical replicates
           alpha.tgt    = 0.80,      # Efficiency of target gene
           alpha.ref    = 0.95,      # Efficiency of reference gene
           sample.sd    = 1,         # Sample standard deviation
           ddcq         = 1) {       # Delta Delta Cq value / effectsize delta

  genes   <- c("tgt", "ref")
  types   <- c("case", "ctrl")
  mu.gene <- c(mu.tgt, mu.ref)
  l <- n.dilutions
  m <- n.replicates
  n <- n.samples

  if (std.curve && (n.dilutions == 0 || n.dilutions == 1)) {
    warning("n.dilutions need to be strictly larger than 1 for meaningful ",
            "standard curves")
    std.curve <- FALSE
  }
  sim.data  <-
    data.frame(sampleName = sprintf("S%03d", rep(rep(1:(2*n), 2), each = m)),
               geneType   = rep(genes, each = 2*n*m),
               sampleType = rep(types, each = n*m, times = 2),
               replicate  = as.character(rep(1:m, 4*n)))
  mus.gene <- rep(mu.gene, each = 2*n*m)
  mus.type <- c(rep(c(ddcq, 0), each = n*m),  rep(0, 2*n*m))
  error    <- rnorm(4*n*m, mean = 0, sd = tech.sd)
  rnd.eff  <- rep(rnorm(2*n, 0, sample.sd), each = m, times = 2)
  effi     <- rep(c(alpha.tgt, alpha.ref), each = 2*n*m)

  sim.data$Cq <- (mus.gene + mus.type)/effi + rnd.eff + error

  if (std.curve) {
    # Simulating standard curve data
    sim.data$copyNumber <- rep(1, 4*n*m)
    sim.data$l2con      <- -log2(sim.data$copyNumber)
    dilution.series     <- 2^rev(-seq_len(l) + 1)

    std.data <-
      data.frame(sampleName = rep("D001", 2*m*l),
                 geneType   = rep(genes, each = l*m),
                 sampleType = rep("Standard", 2*m*l),
                 replicate  = as.character(rep(1:m, 2*l)),
                 copyNumber = rep(dilution.series, each = m, times = 2),
                 l2con      = -log2(rep(dilution.series, each = m, times = 2)))
    error <- rnorm(2*m*l, mean = 0, sd = tech.sd)
    effi <- rep(c(alpha.tgt, alpha.ref), each = l*m)
    std.data$Cq <- 1/effi*(rep(mu.gene, each = l*m) + std.data$l2con) + error

    sim.data <- rbind(sim.data, std.data)
    attr(sim.data, "std.curve") <- TRUE
  } else {
    sim.data$l2con      <- 0
    sim.data$copyNumber <- 1
    attr(sim.data, "std.curve") <- FALSE
  }

  sim.data$geneName <-
    gsub("ref", "GOI", gsub("tgt", "house", as.character(sim.data$geneType)))
  class(sim.data) <- c("data.qPCR", "data.frame")

  return(sim.data)
}


#
# Fitting a qPCRData object
#

qPCRfit <- function(data, weighted, ...) {
  # data is a "data.qPCR" object
  # ... arguments passed to lme

  if (!is.data.qPCR(data)) {
    data <- as.data.qPCR(data)
  }

  # data$RE <- with(data, droplevels(sampleType:sampleName))
  data$RE <- data$sampleName

  if (std.curve(data)) {

    data$std.d <- ifelse(data$sampleType == "Standard", 1, 0)
    data$nstd.d <- 1 - data$std.d
    data$VI <- factor(data$nstd.d)

    if (weighted) {
      weights <- varIdent(form = ~1 | VI)
    } else {
      weights <- NULL
    }

    fit <- lme(Cq ~ -1 + sampleType:geneType + l2con:geneType,
               random = ~1 | RE,
               weights = weights,
               data = data, method = "ML",
               control = lmeControl(returnObject = TRUE), ...)

  } else {

    fit <- lme(Cq ~ -1 + sampleType:geneType,
               random = ~1 | RE,
               data = data, method = "ML", ...)

  }

  return(fit)
}

#
# Delta delta Cq analysis method
#

cHyp <- function(fit, eff.cor) {
  # Function for evaluating the mapping into DDCq
  e <- fixef(fit)
  names(e) <- gsub("sample|gene|Type", "", names(e))
  eff <- e[c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")]

  if (eff.cor) {
    gam <- e[rep(c("tgt:l2con","ref:l2con"), 2)]
  } else {
    gam <- rep(1, length(eff))
  }
  return(sum(gam^-1*eff*c(1, -1, -1, 1)))
}

DcHyp <- function(fit, var.adj) {
  # Compute the gradient of cHyp
  e <- fixef(fit)
  names(e) <- gsub("sample|gene|Type", "", names(e))
  eff <- e[c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")]
  std.data <- any(grepl("Standard|l2con", names(e)))

  if (std.data) {
    M <- cbind(c(-1, 0, 1, 0), c(0, 1, 0, -1))
    gam <- e[rep(c("tgt:l2con","ref:l2con"), 2)]
    D <- c(gam^-1*c(1, -1, -1, 1), (gam^-2*eff) %*% M * (var.adj))
    return(D)
  } else {
    return(c(1, -1, -1, 1))
  }
}

Var.cHyp <- function(fit, var.adj) {
  # Compute variance estimate
  e <- fixef(fit)
  names(e) <- gsub("sample|gene|Type", "", names(e))
  v <- vcov(fit)
  rownames(v) <- colnames(v) <- names(e)

  std.data <- any(grepl("Standard|l2con", names(e)))

  get <- c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref",
           if (std.data) {c("tgt:l2con", "ref:l2con")})

  grad <- DcHyp(fit, var.adj)
  return(as.numeric(t(grad) %*% v[get, get] %*% grad))
}

Var.cHyp.monte.carlo <- function(fit, n = 1e6) {
  # Evalutate the variance by monte carlo integration
  # About fourth decimal place n = 1e5
  e <- fixef(fit)
  names(e) <- gsub("sample|gene|Type", "", names(e))
  v <- vcov(fit)
  rownames(v) <- colnames(v) <- names(e)
  std.data <- any(grepl("Standard|l2con", names(e)))

  rmvn <- rmvnormal(n, e, as.matrix(v))
  colnames(rmvn) <- names(e)

  # Truncate if efficiencies are too low
  low <- abs(rmvn[,c("tgt:l2con","ref:l2con")]) < 0.05
  is.low <- rowSums(low)
  n.low <- sum(is.low)
  if (n.low) {
    warning("Some of the realised efficiencies in monte carlo integration are",
            " very low. Perhaps the std. error of the efficiency is too high.",
            " Removed ", n.low, " values")
    rmvn <- rmvn[!is.low, ]
  }

  if (std.data) {
    gam <- rmvn[, rep(c("tgt:l2con","ref:l2con"), 2)]
  } else {
    gam <- 1
  }

  rmvn <- rmvn[, c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")]
  var.chyp <- var(rowSums(gam^-1*t(t(rmvn)*c(1, -1, -1, 1))))
  return(var.chyp)
}


DDCq <- function(data, eff.cor, var.adj, alpha,
                 var.type = c("deltamethod", "montecarlo"),
                 weighted) {
  # Function to calculate efficiency corrected with or without
  # adjusted variance DDCq values in qPCR experiments
  # data is a qPCR.data object but can also be a fit from lmer
  # var.adj is a boolean value
  # var.type is the method used to estimate the variance of DDCq

  var.type <- match.arg(var.type)

  if (inherits(data, "lmerMod") || inherits(data, "lme")) {
    fit <- data
    std.data <- any(grepl("Standard|l2con", names(fixef(fit))))
  } else if (is.data.qPCR(data)) {
    fit <- qPCRfit(data, weighted = weighted)
    std.data <- std.curve(data)
  } else {
    stop("data should be a data.qPCR dataset or a lmer fit.")
  }

  # Perform t test using above functions
  con <- cHyp(fit, eff.cor)
  if (var.type == "deltamethod") {
    var.con <- Var.cHyp(fit, var.adj)
  }
  if (var.type == "montecarlo") {
    if (!var.adj) warning("If montecarlo")
    var.con <- Var.cHyp.monte.carlo(fit)
  }
  se.con  <- sqrt(var.con)
  t       <- con/se.con

  # Assumes paired samples
  # df <- getME(fit, "q") - 2 - 2 - 2*std.data
  if (inherits(fit, "lme")) df <- unique(summary(fit)$tTable[,"DF"])
  if (inherits(fit, "lmerMod")) df <- getME(fit, "q") - 2 - 2*std.data
  if (length(df) != 1) {
    warning("df is not unique in fit.")
    df <- min(df)
  }

  stopifnot(df > 0)
  p.val   <- 2*(1 - pt(abs(t), df))
  conf.int <- con + c(-1, 1)*qt(1 - alpha/2, df)*se.con
  result  <- c("Estimate" = con, "Std. Error" = se.con,
               "t value" = t, "df" = round(df), "Pr(>|t|)" = p.val,
               "LCL" = conf.int[1], "UCL" = conf.int[2])
  return(result)
}

#
# Wrapper for using DDCq with various methods
#

addCI <- function(x, alpha = 0.05) {
  t <- qt(1 - alpha/2, df = x["df.n"])
  nms <- names(x)
  res <- c(x, x["Estimate"] + c(-1, 1)*t*x["Std. Error"])
  names(res) <- c(nms, "LCL", "UCL")
  return(res)
}

DDCq.test <- function(data,
                      method = c("LMM", "Naive", "Bootstrap", "pBootstrap"),
                      eff.cor = TRUE,
                      var.adj = eff.cor,
                      subset.cols,
                      subset.rows,
                      n.boots = 100,
                      alpha = 0.05,
                      var.type = c("deltamethod", "montecarlo"),
                      weighted = FALSE,
                      ...) {
  # method; character, either "Naive" (i.e a simple t-test) or "LMM"
  # eff.cor; logical, Should the data be efficiency corrected?
  # var.adj; logical, Should the eff. corrected estimate be variace corrected?
  # subset.cols and subset.rows can be used to subset the rows and columns
  method <- match.arg(method)
  var.type <- match.arg(var.type)

  if (method == "LMM") {
    if (eff.cor == FALSE && var.adj == TRUE) {
      eff.cor <- TRUE
      message("eff.cor was coerced to TRUE as var.adj is TRUE.")
    }
    stopifnot(is.logical(eff.cor))
    stopifnot(is.logical(var.adj))
  }

  if (!missing(subset.cols)) {
    data <- data[, subset.cols]
  }

  if (!missing(subset.rows)) {
    data <- data[subset.rows, ]
  }

  if (method == "Naive") {    # Naive method

    data  <- subset(data, l2con == 0 & sampleType != "Standard") # Ignore dil.
    hmean <- aggregate(Cq ~ sampleName + sampleType + geneName + geneType,
                       data = data, FUN = mean)
    wmean <- reshape(hmean, idvar = c("sampleName","sampleType"),
                     timevar = c("geneType"), direction = "wide")
    colnames(wmean) <- gsub("Cq.", "", colnames(wmean))

    t <- t.test(tgt - ref ~ sampleType, var.equal = TRUE, data = wmean,
                conf.level = 1 - alpha)
    est <- -diff(t$estimate)

    result <-
      c(est, est/t$statistic, t$statistic, t$parameter, t$p.value, t$conf.int)
    names(result) <-
      c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)", "LCL", "UCL")
    return(result)

  } else if (method == "LMM") {    # Linear mixed model method

    # Average over replicates
    data <-  aggregate(Cq ~ sampleName + sampleType + geneName + geneType +
                         copyNumber + l2con, data = data, FUN = mean)

    if (eff.cor == FALSE) {
      # Simple DDCq method,  Ignoring dilution data
      data <- subset(data, l2con == 0 & sampleType != "Standard")

      # Note, the use var.adj has no impact
      return(DDCq(as.data.qPCR(data), eff.cor = eff.cor, var.adj = var.adj,
                  alpha = alpha, var.type = var.type, weighted = weighted))

    } else {
      # LMM DDCq method efficiency corr + variance adj. if var.adj == TRUE

      return(DDCq(as.data.qPCR(data), eff.cor = eff.cor, var.adj = var.adj,
                  alpha = alpha, var.type = var.type, weighted = weighted))

    }
  } else if (method == "Bootstrap") {

    return(bootstrapEstimate(data, n.boots, alpha, weighted = weighted))

  } else if (method == "pBootstrap") {

    warning("The parametric boostrap estimate only works with lmerMod objects")
    return(parametricBootstrapEstimate(data, n.boots, alpha, ...))

  } else {

    stop("No usable method found")

  }
}


#
# Bootstrapping method and faclilities
#

Bootstrap.qPCR <- function(data, # A data.qPCR object
                           start.sample   = 3,
                           n.samples      = 5,
                           start.dilution = 3,
                           n.dilutions    = 5,
                           n.resamp       = 100,
                           verbose        = TRUE) {
  # THIS FUNCTION IS NOT USED ANYMORE
  res      <- NULL
  samp.res <- matrix(NA, nrow = 2, ncol = n.resamp)
  dil.res  <- matrix(NA, nrow = 6, ncol = n.samples)
  colnames(dil.res) <- paste("samples =", 1:n.samples + start.sample - 1)
  rownames(dil.res) <- c("mean.est",  "0.05quan.est",  "97.5quan.est",
                         "mean.pval", "0.05quan.pval", "97.5quan.pval")

  # Internal sampling function
  samp <- function(data, sampleType, geneType, j) {
    tmp.data <- data[data$sampleType ==  sampleType &
                       data$geneType   ==  geneType, ]
    nsamp    <- sample(nrow(tmp.data), j, replace = TRUE)
    tmp.data <- tmp.data[nsamp, ]
    tmp.data$sampleName <- paste("H", 1:nrow(tmp.data), sep = "")
    return(tmp.data)
  }

  for (i in seq_len(n.dilutions) + start.dilution - 1) {
    for (j in  seq_len(n.samples) + start.sample - 1) {
      for (k in seq_len(n.resamp)) {

        names.cons <- sort(unique(data$l2con))
        data       <- data[data$l2con %in% names.cons[seq_len(i)], ]

        # Constructing sampled data set
        ntref.data <- samp(data, sampleType = "case", geneType = "tgt", j)
        nrref.data <- samp(data, sampleType = "case", geneType = "ref", j)
        etref.data <- samp(data, sampleType = "ctrl", geneType = "tgt", j)
        erref.data <- samp(data, sampleType = "ctrl", geneType = "ref", j)
        d.data     <- data[data$sampleType == "Standard", ]

        new.data   <- as.data.qPCR(rbind(d.data, ntref.data,
                                         nrref.data, etref.data, erref.data))

        # Efficientcy corrected and variance adjusted method
        tmp.ddcq     <- DDCq.test(new.data, eff.cor = TRUE, var.adj = TRUE)
        samp.res[,k] <- tmp.ddcq[c(1, 5)]
      }

      mean.res     <- rowMeans(samp.res)
      quantile.est <- quantile(samp.res[1,], probs = c(0.25, 0.975))
      quantile.p   <- quantile(samp.res[2,], probs = c(0.25, 0.975))

      dil.res[c(1,4), j - start.sample + 1] <- mean.res
      dil.res[c(2,3), j - start.sample + 1] <- quantile.est
      dil.res[c(5,6), j - start.sample + 1] <- quantile.p

      if (verbose) {
        cat(i, "dilutions", j, "samples bootstrap done.\n")
        flush.console()
      }
    }
    res <- c(res, list(dil.res))
  }
  names(res) <- paste("dilutions =", 1:n.dilutions + start.dilution - 1)
  return(res)
}



bootstrapSample <- function(data) {
  # Function to get one bootstrap sample
  # data should be aggregated

  data$Cq_old <- data$Cq
  std <- data$sampleType == "Standard"

  # Bootstrap the standard dilution curve data
  residualBootstrap <- function(data) {
    lin.fit <- lm(Cq ~ l2con, data = data)
    n <- nrow(data)
    repeat { # make sure
      ind <- sample(seq_len(n), replace = TRUE)
      if (length(unique(ind)) != 1L) break
    }
    new.resid <- resid(lin.fit)[ind]
    return(fitted(lin.fit) + new.resid)
  }

  for (gtype in unique(data$geneType)) {
    if (!is.null(data$geneName)) {
      for (gname in unique(data$geneName)) {
        l <- with(data, std & geneType == gtype & geneName == gname)
        if (sum(l) > 0) {
          data$Cq[l] <- residualBootstrap(data[l, ])
        }
      }
    } else {
      l <- with(data, std & geneType == gtype)
      if (sum(l) > 0) {
        data$Cq[l] <- residualBootstrap(data[l, ])
      }
    }
  }

  # Bootstrap the case/ctrl data
  data$sampleName <- as.character(data$sampleName)
  for (stype in setdiff(unique(data$sampleType), "Standard")) {

      get.tgt <- with(data, geneType == "tgt" & sampleType == stype)
      get.ref <- with(data, geneType == "ref" & sampleType == stype)

      repeat { # make sure
        ind <- sample(seq_len(sum(get.tgt)), replace = TRUE)
        if (length(unique(ind)) != 1L) break
      }

      data[get.tgt, "Cq"] <- data[get.tgt, "Cq_old"][ind]
      data[get.ref, "Cq"] <- data[get.ref, "Cq_old"][ind]

  }
  data$sampleName <- as.factor(data$sampleName)

  data <- subset(data, select = -Cq_old)
  return(as.data.qPCR(data))
}

twoSideP <- function(bdist){
  # From
  #  http://stats.stackexchange.com/questions/83012/
  #   how-to-obtain-p-values-of-coefficients-from-bootstrap-regression
  p1 <- (sum(bdist >= 0) + 1)/(length(bdist) + 1)
  p2 <- (sum(bdist <= 0) + 1)/(length(bdist) + 1)
  p <- 2*min(p1, p2)
  return(p)
}

bootstrapEstimate <- function(data, n.boots, alpha, weighted) {
  # Create bootstrap samples
  bs.samples <- replicate(n.boots, attachWarning(bootstrapSample(data)),
                          simplify = FALSE)

  # Apply DDCq estimate
  res <- sapply(bs.samples, DDCq.test, eff.cor = TRUE, var.adj = FALSE,
                weighted = weighted)

  # Compute the statistics
  ddcq <- res["Estimate", ]

  ans <- c("Estimate" = mean(ddcq), "Std. Error" = sd(ddcq),
           "t value" = NA, "df" = NA, "Pr(>|t|)" = twoSideP(ddcq),
           "LCL" = quantile(ddcq, alpha/2), "UCL" = quantile(ddcq, 1 - alpha/2))

  attr(ans, "extra") <- res
  # attr(ans, "warnings") <- lapply(bs.samples, function(x) attr(x, "warnings"))

  return(ans)
}


parametricBootstrapEstimate <- function(data, n.boots, alpha, weighted) {

  fit <- qPCRfit(data, weighted = weighted)
  ddcq <- function(x) {
    catchBootstrapWarning(DDCq(x, var.adj = FALSE, alpha = alpha,
                               eff.cor = TRUE, var.adj = FALSE,
                               weighted = weighted)["Estimate"])
  }

  r <- bootMer(fit, ddcq, nsim = n.boots)
  t <- r$t[!is.na(r$t)] # Omit NAs

  ans <- c("Estimate" = mean(t), "Std. Error" = sd(t),
           "t value" = NA, "df" = NA, "Pr(>|t|)" = twoSideP(t),
           "LCL" = quantile(t, alpha/2), "UCL" = quantile(t, 1 - alpha/2))

  attr(ans, "extra") <- r
  return(ans)
}



#
# Function to coerce data.frame into data.qPCR object
#

as.data.qPCR <- function(data) {
  if (!suppressMessages(is.data.qPCR(data))) {
    class(data) <- c("data.qPCR", "data.frame")

    if (!"replicate" %in% names(data)) {
      data$replicate <- 1
    }
    if ("Standard" %in% data$sampleType || length(unique(data$l2con)) != 1L) {
      attr(data, "std.curve") <- TRUE
    } else {
      attr(data, "std.curve") <- FALSE
    }
    if (!is.data.qPCR(data)) {
      stop("Could not coerce data.frame to data.qPCR object")
    }
    data$geneType   <- as.factor(as.character(data$geneType))
    data$sampleType <- as.factor(as.character(data$sampleType))
    data$sampleName <- as.factor(as.character(data$sampleName))
    if ("geneName" %in% names(data)) {
      data$geneName <- as.factor(as.character(data$geneName))
    }
  }
  return(data)
}



################################################################################
# Auxilary functions                                                           #
################################################################################

#
# Function to get and set "std.curve" attribute
#

std.curve <- function(data) {
  if (!is.data.qPCR(data))
    stop("data is not a data.qPCR object.")
  attr(data, "std.curve")
}

`std.curve<-` <- function(value, data) {
  if (!is.data.qPCR(data))
    stop("data is not a data.qPCR object.")
  data <- structure(data, std.curve = value)
}

#
# Function for regularising power estimates (convex hull)
#

cenvelope <- function(x) {
  # For regularising the power calculations. x is a 2 by n matrix.
  x <- x[chull(x), ]
  if (x[1,1] != min(x[, 1])) {
    x <- x[which(x[, 1] == min(x[, 1])):length(x[, 1]), ]
  }
  return(x)
}

#
# Function for checking whether an object is of data.qPCR class
#

is.data.qPCR <- function(object) {
  arg <- deparse(substitute(object))
  if (!is.data.frame(object)) {
    message(arg, "is not a data.frame.\n")
    return(FALSE)
  }
  if (!identical(c("data.qPCR", "data.frame"), class(object))) {
    message("Class of ", arg, " is not correctly defined.\n")
    return(FALSE)
  }
  if (is.null(attr(object, "std.curve"))) {
    message(arg, " has no logical std.curve attribute.\n")
    return(FALSE)
  }
  cols <- c("sampleName", "geneType", "sampleType",
            "replicate", "Cq", "copyNumber", "l2con")
  if (attr(object, "std.curve")) {
    if (!all(cols %in% colnames(object))) {
      message("Needed columns of ", arg, " are not present.\n",
              "The missing column(s) are:\n",
              cols[!(cols %in% colnames(object))])
      return(FALSE)
    }
  } else {
    if (!all(cols[-6:-7] %in% colnames(object))) {
      message("Needed columns of ", arg, " are not present.\n",
              "The missing column(s) are:\n",
              (cols[-6:-7])[!(cols[-6:-7] %in% colnames(object))])
      return(FALSE)
    }
  }
  return(TRUE)
}

#
# sort qPCR data
#

sort.data.qPCR <- function(data) {
  data <- data[order(data$sampleType), ]
  data <- data[order(data$geneType), ]
  if (std.curve(data)) {
    data <- data[order(data$sampleType == "Standard"), ]
  }
  return(data)
}


#
# Formatting p-values function
#

sn <- function(x, digits) {
  ord <- floor(log10(abs(x)))
  y   <- x/10^ord
  if (!missing(digits)) {
    y <- format(y, digits = digits)
  }
  y <- paste(format(y), "\\cdot 10^{", ord, "}", sep = "")
  if (any(ord %in% c(-Inf, 0))) {
    y[ord %in% c(-Inf, 0)] <- as.character(x)[ord %in% c(-Inf, 0)]
  }
  return(paste("$", y, "$", sep = ""))
}

#
# Add to .RData file
#

resave <- function(..., list = character(), file) {
  # Resave function from http://github.com/AEBilgrau/Bmisc
  if (!file.exists(file)) { # If file does not exists resave functions as save
    save(..., list = list, file = file)
  }
  previous  <- load(file)  # Returns the loaded object names
  var.names <- c(list, as.character(substitute(list(...)))[-1L])
  for (var in var.names) {
    assign(var, get(var, envir = parent.frame()))
  }
  save(list = unique(c(previous, var.names)), file = file)
}


attachWarning <- function(expr) {
  # Function to handle  warnings (attach as attribute)
  w.handler <- function(w) { # warning handler
    W <<- w
    invokeRestart("muffleWarning")
  }
  W <- NULL
  value <- withCallingHandlers(tryCatch(expr), warning = w.handler)
  attr(value, "warnings") <- W
  return(value)
}


catchBootstrapWarning <- function(expr) {
  warningHandler <- function(w) {
    out <- NA
    attr(out, "warning") <- w
    return(out)
  }
  return(tryCatch(expr, warning = warningHandler))
}

catchBootstrapWarning(log(-1))

export <- # Function names to export to parallel clusters
  c("SimqPCRData", "qPCRfit", "is.data.qPCR", "cHyp", "DcHyp", "Var.cHyp",
    "std.curve", "DDCq.test", "DDCq", "as.data.qPCR", "bootstrapEstimate",
    "parametricBootstrapEstimate", "bootstrapSample", "twoSideP",
    "catchBootstrapWarning", "attachWarning", "Var.cHyp.monte.carlo")





