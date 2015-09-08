
################################################################################
# Analysis of Yuan data                                                        #
# Written by:                                                                  #
#   Anders Ellern Bilgrau                                                      #
################################################################################


#
# Plot Yuan et al. (2008) data
#

lime.colours <- c("#C4B582", "#65654D", "#B3938C")

mypanel <- function(x, y, ...) {

  panel.lmlineq(x, y, adj = c(1,0), lty = 2, col.line = "darkgrey", digits = 2,
                at = 0.5, rot = TRUE, cex = 0.7, style = 2, pos = 1,
                varNames = alist(y = C[q], x = k), offset = 1)
  panel.xyplot(x, y, ...)
  aa <- lm(y ~ x)
  panel.text(2, 28, labels = sprintf("se = %.3f", sqrt(vcov(aa)[2,2])),
             cex = 0.7)
}

fig3 <- xyplot(Cq ~ l2con | sampleType:geneName, data = yuan, panel = mypanel,
               groups = yuan$geneName,
               col = lime.colours, pch  = seq_along(lime.colours),
               xlab = "Dilution step", #as.expression(bquote(-log[2]*N["0,i,j,k"])),
               ylab = expression(C[q]),
               main = "")

setEPS()
postscript("../output/fig3.eps", width = 1.5*7, height = 0.5*7, fonts = "serif")
  trellis.par.set(strip.background = list(col = "lightgrey"))
  print(fig3, position = c(0, 0, 1, 1))
dev.off()


#
# Analysis
#

DDCq.test2 <- function(data, eff.cor = TRUE, var.adj = TRUE, alpha = 0.05) {
  data$re <-
    factor(with(data, paste(sampleType, geneType, copyNumber, sep = ":")))

  fit <- lme(Cq ~ -1 + sampleType:geneType + geneType:sampleType:l2con,
             random = ~1 | re, data = data, method = "ML")

#   data <- aggregate(Cq ~ ., FUN = mean, data = data)
#   fit <- lm(Cq ~ -1 + sampleType:geneType + geneType:sampleType:l2con,
#             data = data)
#   fixef.lm <- function(x) {coef(x)}

  # Function for evaluating the mapping into DDCq and it variance
  e <- fixef(fit)
  v <- vcov(fit)
  names(e) <- gsub("sample|gene|Type", "", names(e))
  names(e) <- gsub("(tgt|ref):(case|ctrl)", "\\2:\\1", names(e))
  rownames(v) <- colnames(v) <- names(e)

  # Is standard curve data present?
  std.e <- grepl("Standard|l2con", names(e))
  std.data <- any(std.e)

  # Order and prep names
  sign <- c(1, -1, -1, 1)
  mu.nms <- c("case:tgt", "case:ref", "ctrl:tgt", "ctrl:ref")
  gam.nms <- paste0(mu.nms, ":l2con")
  if (!all(gam.nms %in% names(e))) {
    gam.nms <- gsub("case:|ctrl:", "", gam.nms)
  }
  stopifnot(all(gam.nms %in% names(e)))

  e <- e[c(mu.nms, gam.nms)]
  v <- v[c(mu.nms, gam.nms), c(mu.nms, gam.nms)]

  # Get mus in correct order
  mu <- e[mu.nms]

  # Get gammas in correct order
  gam <- e[gam.nms]
  if (!std.data || !eff.cor) {
    gam[] <- 1
  }

  # compute ddcq
  con <- sum(mu*sign*gam^(-1))

  # Compute gradient in estimate
  D <- c(gam*sign, -sign*gam^(-2)*mu*var.adj)
  var.con <- c(tcrossprod(D, v) %*% D)

  # Generate results
  se.con <- sqrt(var.con)
  t.stat <- con/se.con
  # df <- fit$df.residual
  df <- unique(summary(fit)$tTable[,"DF"])
  p.val   <- 2*(1 - pt(abs(t.stat), df))
  conf.int <- con + c(-1, 1)*qt(1 - alpha/2, df)*se.con
  result  <- c("Estimate" = con, "Std. Error" = se.con,
               "t value" = t.stat, "df" = round(df), "Pr(>|t|)" = p.val,
               "LCL" = conf.int[1], "UCL" = conf.int[2])

  return(result)
}


#
# Fit model
#

g1 <- c("MT7", "Tublin")
g2 <- c("MT7", "UBQ")
ds1 <- subset(yuan, geneName %in% g1)
ds2 <- subset(yuan, geneName %in% g2)

# Assuming different slope for each i, j.
toTeX <- rbind(EC    = DDCq.test2(ds1, eff.cor = TRUE, var.adj = FALSE),
               ECVA1 = DDCq.test2(ds1, eff.cor = TRUE, var.adj = TRUE),
               EC    = DDCq.test2(ds2, eff.cor = TRUE, var.adj = FALSE),
               ECVA1 = DDCq.test2(ds2, eff.cor = TRUE, var.adj = TRUE))



toTeX      <- signif(toTeX, 4)
toTeX[, 5] <- sn(toTeX[, 5])
colnames(toTeX) <- gsub("Pr(>|t|)", "$p$-value", colnames(toTeX), fixed = TRUE)
colnames(toTeX) <- gsub("t ", "$t$-", colnames(toTeX), fixed = TRUE)

rownames(toTeX) <- gsub("LMEM", "LMM", rownames(toTeX))
rownames(toTeX) <- gsub("t.", "$t$-", rownames(toTeX), fixed = TRUE)
rownames(toTeX) <- gsub("ECVA", "EC\\\\&VA", rownames(toTeX))

grps <- sapply(list(g1, g2), function(x) paste(x[1], "vs", x[2]))

caption.txt <- "\\citet{Yuan2008} data: Method comparison for estimating the
  $\\ddcq$-value. EC denotes use of the plug-in estimator disregarding
  uncertainty in the AE.
  EC\\&VA1 denotes that the efficiency correction was variance adjusted using
  the delta method."
w <- latex(toTeX,
           file    = "../output/Table3.tex",
           title   = "",
           label   = "table:yuan",
           caption = sprintf(caption.txt, n.boots, "\\%"),
           caption.loc = "top",
           rgroup  = grps,
           center  = "center",
           numeric.dollar = TRUE,
           keep.tex = TRUE,
           size = "small")
