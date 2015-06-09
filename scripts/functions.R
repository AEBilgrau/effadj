
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
           n.samples    = 8,         # Number of samples
           n.replicates = 3,         # Number of replicates
           n.dilutions  = 5,         # Number of dilutions
           tech.sd      = 1/2,       # Variance of technical replicates
           alpha.tgt    = 0.9,       # Efficiency of target gene
           alpha.ref    = 0.9,       # Efficiency of reference gene
           sample.sd    = 1,         # Sample standard deviation
           ddcq         = 1) {       # Delta Delta Cq value / effectsize delta

  genes   <- c("tgt", "ref")
  types   <- c("case", "ctrl")
  mu.gene <- c(mu.tgt, mu.ref)
  l <- n.dilutions
  m <- n.replicates
  n <- n.samples

  if (n.dilutions == 0 | n.dilutions == 1) {
    std.curve <- FALSE
  }
  sim.data  <-
    data.frame(sampleName = sprintf("S%03d", rep(rep(1:n, 4), each = m)),
               geneType   = rep(genes, each = 2*n*m),
               sampleType = rep(types, each = n*m, times = 2),
               replicate  = as.character(rep(1:m, 4*n)))
  mus.gene <- rep(mu.gene, each = 2*n*m)
  mus.type <- c(rep(c(ddcq, 0), each = n*m),  rep(0, 2*n*m))
  error    <- rnorm(4*n*m, mean = 0, sd = tech.sd)
  rnd.eff  <- rep(rnorm(2*n, 0, sample.sd), each = m, times = 2)

  sim.data$Cq <-
    rep(c(alpha.tgt, alpha.ref), each = 2*n*m)^-1*
    (mus.gene + mus.type) + rnd.eff + error

  if (std.curve) {
    # Simulating standard curve data
    sim.data$copyNumber <- rep(1, 4*n*m)
    sim.data$l2con      <- -log2(sim.data$copyNumber)
    dilution.series     <- 2^((1-l):0)

    std.data <-
      data.frame(sampleName = rep("D001", 2*m*l),
                 geneType   = rep(genes, each = l*m),
                 sampleType = rep("Standard", 2*m*l),
                 replicate  = as.character(rep(1:m, 2*l)),
                 copyNumber = rep(dilution.series, each = m, times = 2),
                 l2con      = -log2(rep(dilution.series, each = m, times = 2)))
    error       <- rnorm(2*m*l, mean = 0, sd = tech.sd)
    std.data$Cq <-
      rep(c(alpha.tgt, alpha.ref), each = l*m)^-1*
      (rep(mu.gene, each = l*m) + std.data$l2con) + error

    sim.data <- rbind(sim.data, std.data)
    attr(sim.data, "std.curve") <- TRUE
  } else {
    sim.data$l2con      <- 0
    sim.data$copyNumber <- 1
    attr(sim.data, "std.curve") <- FALSE
  }
  class(sim.data) <- c("data.qPCR", "data.frame")
  return(sim.data)
}


#
# Fitting a qPCRData object
#

qPCRfit <- function(data, ...) {
  # data is a "data.qPCR" object
  # ... arguments passed to lmer

  if (!is.data.qPCR(data)) {
    arg <- deparse(substitute(data))
    stop(paste(arg, "is not of data.qPCR class."))
  }

  if (std.curve(data)) {

    fit <- lmer(Cq ~ -1 + sampleType:geneType + l2con:geneType +
                  (1 | sampleType/sampleName),
                data = data, REML = FALSE, ...)

  } else {

    fit <- lmer(Cq ~ -1 + sampleType:geneType + (1 | sampleType/sampleName),
                data = data, REML = FALSE, ...)

  }
  return(fit)
}

#
# Delta delta Cq analysis method
#

DDCq <- function (data, var.adj) {
  # Function to calculate efficiency corrected with or without
  # adjusted variance DDCq values in qPCR experiments
  # data is a qPCR.data object
  # var.adj is a boolean value

  # Internal functions:
  cHyp <- function (fit) {
    # Function for evaluating the mapping into DDCq
    e <- fixef(fit)
    names(e) <- gsub("sample|gene|Type", "", names(e))
    eff <- e[c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")]

    if (std.curve(data)) {
      gam <- e[rep(c("tgt:l2con","ref:l2con"), 2)]
    } else {
      gam <- rep(1, length(eff))
    }
    return(sum(eff*gam^-1*c(1, -1, -1, 1)))
  }

  DcHyp <- function (fit, var.adj) {
    # Compute the gradient of cHyp
    e <- fixef(fit)
    names(e) <- gsub("sample|gene|Type", "", names(e))
    eff <- e[c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")]
    if (std.curve(data)) {
      gam <- e[rep(c("tgt:l2con","ref:l2con"), 2)]
      D <- c(gam^-1*c(1, -1, -1, 1),
             (gam^-2*eff) %*% cbind(c(-1, 0, 1, 0), c(0, 1, 0, -1))*(var.adj))
      return(D)
    } else {
      return(c(1, -1, -1, 1))
    }
  }

  Var.cHyp <- function (fit, var.adj = var.adj) {
    # Compute variance estimate
    e <- fixef(fit)
    names(e) <- gsub("sample|gene|Type", "", names(e))
    v <- vcov(fit)
    rownames(v) <- colnames(v) <- names(e)
    get <- c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref",
             if (std.curve(data)) {c("tgt:l2con", "ref:l2con")})
    grad <- DcHyp(fit, var.adj)
    return(as.numeric(t(grad)%*%v[get, get]%*%grad))
  }

  # Perform t test using above functions
  fit     <- qPCRfit(data)
  con     <- cHyp(fit)
  var.con <- Var.cHyp(fit, var.adj)
  t       <- con/sqrt(var.con)
  df      <- getME(fit, "q") - 2 - 2 - 2*std.curve(data)
  p.val   <- 2*(1 - pt(abs(t), df))
  result  <- c("Estimate" = con, "Std. Error" = sqrt(var.con),
               "t value" = t, "df" = round(df), "Pr(>|t|)" = p.val)
  return(result)
}

#
# Wrapper for using DDCq with various methods
#

DDCq.test <- function (data,
                       method = c("LMM", "Naive", "Bootstrap"),
                       eff.cor,
                       var.adj,
                       subset.cols,
                       subset.rows,
                       n.boots = 100) {
  # method; character, either "Naive" (i.e a simple t-test) or "LMM"
  # eff.cor; logical, Should the data be efficiency corrected?
  # var.adj; logical, Should the eff. corrected estimate be variace corrected?
  # subset.cols and subset.rows can be used to subset the rows and columns
  method <- match.arg(method)

  if (!missing(var.adj) & !missing(var.adj)) {
    if (var.adj == TRUE & eff.cor == FALSE) {
      eff.cor <- TRUE
      message("eff.cor was coerced to TRUE as var.adj is TRUE.")
    }
  }
  if (!missing(subset.cols)) {
    data <- data[, subset.cols]
  }
  if (!missing(subset.rows)) {
    data <- data[subset.rows, ]
  }

  if (method == "Naive") {    # Naive method

    data  <- data[data$sampleType != "Standard", ]  # Ignoring dilution data
    hmean <- aggregate(Cq ~ sampleName + geneType + sampleType,
                       data = data, FUN = mean)
    wmean <- reshape(hmean, idvar = c("sampleName","sampleType"),
                     timevar = c("geneType"), direction = "wide")
    colnames(wmean) <- gsub("Cq.", "", colnames(wmean))

    t <- t.test(tgt - ref ~ sampleType, var.equal = TRUE, data = wmean)
    est <- -diff(t$estimate)
    result <- c(est, est/t$statistic, t$statistic, t$parameter, t$p.value)
    names(result) <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
    return(result)

  } else if (method == "LMM") {    # Linear mixed model method

    if (eff.cor == FALSE) {
      # Simple DDCq method
      data <- data[data$sampleType != "Standard", ]  # Ignoring dilution data
      data <- aggregate(Cq ~ sampleName + geneType + sampleType +
                          l2con + copyNumber,
                        data = data, FUN = mean)
      # Note, the use of var.adj has no impact here!
      return(DDCq(as.data.qPCR(data)))
    } else {
      # LMM DDCq method. efficiency corrected and
      # variance adjusted if var.adj == TRUE
      data <- aggregate(Cq ~ sampleName + geneType + sampleType +
                          l2con + copyNumber,
                        data = data, FUN = mean)  # Mean over replicates
      return(DDCq(as.data.qPCR(data), var.adj = var.adj))
    }
  } else if (method == "Bootstrap") {    # Booststrap

    return(bootstrapEstimate(data, n.boots))

  } else {
    stop("No usable method found")
  }
}


#
# Power simulation of standard curves
#

PowerSim <-
  function (n.sims         = 400,
            start.sample   = 5,
            n.samples      = 2,
            start.dilution = 3,
            n.dilutions    = 2,
            alpha.lvl      = 0.05,
            ddcq           = 1.0,
            std.curve      = TRUE,
            method         = "LMM",
            eff.cor        = TRUE,
            var.adj        = TRUE,
            ... ) {   # ... passed to SimqPCRData function
    st <- proc.time()
    if (std.curve == FALSE) {
      start.dilution <- 1
      n.dilutions    <- 1
    }

    pow.res           <- matrix(0, ncol = n.samples, nrow = n.dilutions)
    colnames(pow.res) <- paste("samples =",   1:n.samples   + start.sample  -1)
    rownames(pow.res) <- paste("dilutions =", 1:n.dilutions + start.dilution-1)

    for (k in seq_len(n.dilutions)) {
      tests <- matrix(0, nrow = n.sims, ncol = n.samples)

      for (i in seq_len(n.samples)) {
        cat("Computing power with", start.dilution + k - 1,
            "dilutions and", start.sample + i - 1, "samples.\n")
        flush.console()
        for (j in seq_len(n.sims)) {
          sim.data <- SimqPCRData(std.curve   = std.curve,
                                  n.samples   = start.sample   + i - 1,
                                  n.dilutions = start.dilution + k - 1,
                                  ddcq        = ddcq,
                                  ... )
          sim.ddcq <- DDCq.test(sim.data, method = method,
                                eff.cor = eff.cor, var.adj = var.adj)

          tests[j,i] <- ifelse(sim.ddcq[5] < alpha.lvl, 1, 0)
        }
      }
      pow.res[k, ] <- colMeans(tests)
    }
    cat("Simulation finished in", (proc.time()[3] - st[3])%/%60,
        "minutes.\n"); flush.console()
    return(pow.res)
  }


#
# Bootstrapping method
#

Bootstrap.qPCR <- function(data,               # A data.qPCR object
                           start.sample   = 3,
                           n.samples      = 5,
                           start.dilution = 3,
                           n.dilutions    = 5,
                           n.resamp       = 100,
                           verbose        = TRUE) {
  res      <- NULL
  samp.res <- matrix(NA, nrow = 2, ncol = n.resamp)
  dil.res  <- matrix(NA, nrow = 6, ncol = n.samples)
  colnames(dil.res) <- paste("samples =", 1:n.samples + start.sample - 1)
  rownames(dil.res) <- c("mean.est",  "0.05quan.est",  "97.5quan.est",
                         "mean.pval", "0.05quan.pval", "97.5quan.pval")

  # Internal sampling function
  samp <- function (data, sampleType, geneType, j) {
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
    return(fitted(lin.fit) + sample(resid(lin.fit), replace = TRUE))
  }

  # if (is.null(unique(data$geneName))) warning("no geneName col present in data")
  for (gname in unique(data$geneName)) {
    for (gtype in unique(data$geneType)) {
      if (is.null(gname)) {
        get <- with(data, std & geneType == gtype)
      } else {
        get <- with(data, std & geneType == gtype & geneName == gname)
      }
      if (sum(get) > 0) {
        data$Cq[get] <- residualBootstrap(data[get, ])
      }
    }
  }

  # Bootstrap the case/ctrl data
  data$sampleName <- as.character(data$sampleName)
  for (stype in setdiff(unique(data$sampleType), "Standard")) {
    get <- !std & data$sampleType == stype
    nms <- sample(unique(data$sampleName[get]), replace = TRUE)
    nms <- paste0(nms, "_BSS", seq_along(nms))

    ind <- numeric()
    new.sampleName <- character()
    for (nm in nms) {
      i <- which(get & data$sampleName == gsub("_BSS[0-9]+$", "", nm))
      ind <- c(ind, i)
      new.sampleName <- c(new.sampleName, rep(nm, length(i)))
    }
    data[get, ] <- data[ind, ]
    data$sampleName[get] <- new.sampleName
  }
  data$sampleName <- as.factor(data$sampleName)

  return(data)
}

twoSideP <- function(bdist){
  # From
  #  http://stats.stackexchange.com/questions/83012/
  #   how-to-obtain-p-values-of-coefficients-from-bootstrap-regression
  p1 <- sum(bdist > 0)/length(bdist)
  p2 <- sum(bdist < 0)/length(bdist)
  p <- 2*min(p1, p2)
  return(p)
}

bootstrapEstimate <- function(data, n.boots) {
  # Aggregate data
  if (is.null(data$geneName)) {
    data <- aggregate(Cq ~ sampleName + geneType + sampleType +
                        geneType + l2con + copyNumber,
                      data = data, FUN = mean)
  } else {
    data <- aggregate(Cq ~ sampleName + geneType + sampleType +
                        geneType + l2con + copyNumber + geneName,
                      data = data, FUN = mean)
  }

  # Create bootstrap samples
  bs.samples <- replicate(n.boots, bootstrapSample(data), simplify = FALSE)

  # Apply DDCq estimate
  res <- sapply(bs.samples, DDCq.test, eff.cor = TRUE, var.adj = FALSE)

  # Compute the statistics
  ddcq <- res["Estimate", ]

  ans <- c("Estimate" = mean(ddcq), "Std. Error" = sd(ddcq),
           "t value" = NA, "df" = NA, "Pr(>|t|)" =   twoSideP(ddcq))
  return(ans)
}


#
# Function to coerce data.frame into data.qPCR object
#

as.data.qPCR <- function (data) {
  if (!suppressMessages(is.data.qPCR(data))) {
    class(data) <- c("data.qPCR", "data.frame")

    if (!"replicate" %in% names(data)) {
      data$replicate <- 1
    }
    if ("Standard" %in% data$sampleType) {
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
      data$geneName   <- as.factor(as.character(data$nameName))
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

std.curve <- function (data) {
  if (!is.data.qPCR(data))
    stop("data is not a data.qPCR object.")
  attr(data, "std.curve")
}

"std.curve<-" <- function (value, data) {
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
