# Reproduces Section 5 of the TerravAI paper (multi-centre runoff-plot data).
# The data are not distributed with the package. Usage:
#   TERRAVAI_PLOTS=plots.csv TERRAVAI_OUT=outdir Rscript field_analysis.R
# plots.csv columns: id, centre, soil_loss, runoff, rain, slope, cultivation, lulc, crop, measure, N, P, K
inp <- Sys.getenv("TERRAVAI_PLOTS", "plots.csv"); outd <- Sys.getenv("TERRAVAI_OUT", tempdir())
dir.create(file.path(outd, "out"), showWarnings = FALSE)
suppressMessages(library(terravai))
set.seed(1)
d <- read.csv(inp, stringsAsFactors = FALSE)
d <- d[!is.na(d$soil_loss), ]
d$cult <- ifelse(grepl("along", tolower(d$cultivation)), "along",
                 ifelse(grepl("cross|aross", tolower(d$cultivation)), "across", "unknown"))
lu <- tolower(paste(d$lulc, d$crop))
d$lu <- ifelse(grepl("forest|agroforest|grass|alley|silvi", lu), "tree_grass",
          ifelse(grepl("fallow|degraded", lu), "fallow_degraded",
          ifelse(grepl("hort|orange|turmeric|tumeric|cassava|tea|pepper", lu), "horticulture", "field_crop")))
d$ctl <- grepl("control|cont-|without|along slope|^-$|cultivated fallow", tolower(d$measure)) | is.na(d$measure)
for (v in c("runoff", "rain", "slope")) { d[[paste0(v, "_na")]] <- is.na(d[[v]]); d[[paste0(v, "_f")]] <- ifelse(is.na(d[[v]]), -1, d[[v]]) }
d$cult <- factor(d$cult); d$lu <- factor(d$lu)
d$x <- seq_len(nrow(d)); d$y <- 0
centres <- sort(unique(d$centre))
cat("rows with soil loss:", nrow(d), "\n"); print(table(d$centre)); print(table(d$lu))
res <- list(n = nrow(d), by_centre = table(d$centre))

# Piecewise empirical form refitted on training centres
pw_fit <- function(tr, nut) {
  brk <- c(0, 5, 10, 20, Inf); co <- NULL
  for (k in 1:4) {
    s <- tr[tr$soil_loss >= brk[k] & tr$soil_loss < brk[k + 1] & !is.na(tr[[nut]]), ]
    co <- rbind(co, if (nrow(s) >= 3) coef(lm(s[[nut]] ~ s$soil_loss)) else c(NA, NA))
  }
  list(brk = brk, co = co)
}
pw_pred <- function(f, A) { k <- findInterval(A, f$brk, rightmost.closed = TRUE); k[k > 4] <- 4
  b <- f$co[k, , drop = FALSE]; pmax(b[, 1] + b[, 2] * A, 0) }

form <- function(nut) as.formula(paste(nut, "~ soil_loss + runoff_f + rain_f + slope_f + runoff_na + rain_na + slope_na + cult + lu"))
nuts <- c("N", "P", "K")
M <- 500
out <- list()
for (nut in nuts) {
  dd <- d[!is.na(d[[nut]]), ]
  pub <- dm_nutrient_empirical(dd$soil_loss, nut, "soil")
  pub <- if (inherits(pub, "dm_predictive")) pub$mean else pub
  pw <- rep(NA, nrow(dd)); D <- matrix(NA, nrow(dd), M)
  for (cc in centres) {
    te <- dd$centre == cc; tr <- dd[!te, ]
    if (!any(te)) next
    f <- pw_fit(tr, nut); pw[te] <- pw_pred(f, dd$soil_loss[te])
    q <- suppressWarnings(dm_fit_qrf(form(nut), tr, dd[te, ], copula = FALSE, n_draws = M, seed = 1,
                                     num.trees = 500))
    D[te, ] <- q$draws
  }
  qrf <- dm_draws(D)
  out[[nut]] <- list(obs = dd[[nut]], published = pub, refit = pw, qrf = qrf, centre = dd$centre, soil = dd$soil_loss)
  cat("\n=====", nut, "loss (kg/ha), n =", nrow(dd), "\n")
  cat("-- published SoilRunNutri equations (in-sample, optimistic)\n"); print(dm_validate(dd[[nut]], pub))
  cat("-- same piecewise form, refitted leaving one centre out\n"); print(dm_validate(dd[[nut]], pw))
  cat("-- quantile regression forest, leaving one centre out\n"); print(dm_validate(dd[[nut]], qrf))
}
saveRDS(list(d = d, out = out), file.path(outd, "out", "nutrients.rds"))
suppressMessages(library(terravai))
x <- readRDS(file.path(outd, "out", "nutrients.rds")); d <- x$d
centres <- sort(unique(d$centre))
# Hierarchical log-ratio model: log Y = log A + mu_c + e,  mu_c ~ N(m, tau2), e ~ N(0, s2)
fit_prior <- function(tr, nut) {
  ok <- tr$soil_loss > 0.5 & tr[[nut]] > 0 & !is.na(tr[[nut]])
  y <- log(tr[[nut]][ok]); la <- log(tr$soil_loss[ok]); cc <- factor(tr$centre[ok])
  fe <- lm(y ~ la + cc)                       # within-centre slope b
  b <- unname(coef(fe)["la"])
  r <- y - b * la
  mc <- tapply(r, cc, mean); nc <- tapply(r, cc, length)
  s2 <- sum(residuals(fe)^2) / fe$df.residual
  tau2 <- max(var(mc) - mean(s2 / nc), 0.01)
  list(m = mean(mc), tau2 = tau2, s2 = s2, b = b)
}
posterior <- function(pr, rloc) {   # conjugate update of mu from local log ratios
  k <- length(rloc); if (k == 0) return(c(pr$m, pr$tau2))
  v <- 1 / (1 / pr$tau2 + k / pr$s2); c(v * (pr$m / pr$tau2 + sum(rloc) / pr$s2), v)
}
pred_draws <- function(A, post, s2, M = 1000) {
  mu <- rnorm(M, post[1], sqrt(post[2]))                      # common to all plots
  e <- matrix(rnorm(length(A) * M, 0, sqrt(s2)), length(A), M)
  exp(log(pmax(A, 0.05)) + e + matrix(mu, length(A), M, byrow = TRUE))
}
nse <- function(o, p) 1 - sum((p - o)^2) / sum((o - mean(o))^2)
set.seed(11)
ks <- c(0, 1, 3, 5, 10); R <- 40
tab <- NULL; covtab <- NULL; pooled <- list()
for (nut in c("N", "P", "K")) {
  dd <- d[!is.na(d[[nut]]) & d$soil_loss > 0.5 & d[[nut]] > 0, ]
  for (k in ks) {
    o_all <- p_all <- c(); cov50 <- cov90 <- c()
    for (cc in centres) {
      te <- dd[dd$centre == cc, ]; if (nrow(te) < k + 5) next
      pr <- fit_prior(dd[dd$centre != cc, ], nut)
      for (rep in seq_len(if (k == 0) 1 else R)) {
        cal <- if (k == 0) integer(0) else sample(nrow(te), k)
        ev <- if (k == 0) te else te[-cal, ]
        post <- posterior(pr, log(te[[nut]][cal]) - pr$b * log(te$soil_loss[cal]))
        D <- pred_draws(ev$soil_loss^pr$b, post, pr$s2, 400)
        pm <- exp(pr$b * log(pmax(ev$soil_loss, 0.05)) + post[1])   # predictive median
        o_all <- c(o_all, ev[[nut]]); p_all <- c(p_all, pm)
        lo50 <- apply(D, 1, quantile, .25); hi50 <- apply(D, 1, quantile, .75)
        lo90 <- apply(D, 1, quantile, .05); hi90 <- apply(D, 1, quantile, .95)
        cov50 <- c(cov50, ev[[nut]] >= lo50 & ev[[nut]] <= hi50)
        cov90 <- c(cov90, ev[[nut]] >= lo90 & ev[[nut]] <= hi90)
      }
    }
    lo <- log(o_all); lp <- log(p_all)
    tab <- rbind(tab, data.frame(nutrient = nut, k = k, NSE = nse(o_all, p_all), NSE_log = nse(lo, lp),
                                 PBIAS = 100 * sum(p_all - o_all) / sum(o_all),
                                 cov50 = mean(cov50), cov90 = mean(cov90)))
  }
}
print(tab, digits = 3)
for (nut in c('N','P','K')) { pr <- fit_prior(d[!is.na(d[[nut]]) & d$soil_loss > 0.5 & d[[nut]] > 0, ], nut); cat(nut, 'b =', round(pr$b,2), ' tau =', round(sqrt(pr$tau2),2), ' sigma_w =', round(sqrt(pr$s2),2), ' share of variance between centres =', round(pr$tau2/(pr$tau2+pr$s2),2), '\n') }
saveRDS(tab, file.path(outd, "out", "calibration_curve.rds"))
suppressMessages(library(terravai))
x <- readRDS(file.path(outd, "out", "nutrients.rds")); d <- x$d
nut <- "N"; Tn <- 10; r <- 3; k <- 5; R <- 40
dd <- d[!is.na(d$N) & d$soil_loss > 0.5 & d$N > 0, ]
loss_log <- dm_loss(log(Tn), "above", ratio = r)
regret <- function(act, y) mean(ifelse(act & y <= Tn, 1, 0) + ifelse(!act & y > Tn, r, 0))
decide <- function(mu, v, y_known = NULL, idx = NULL) {
  p <- 1 - pnorm((log(Tn) - mu) / sqrt(v)); act <- p > 1 / (1 + r)
  if (length(idx)) act[idx] <- y_known > Tn
  act
}
set.seed(21); res <- NULL
for (cc in names(which(table(dd$centre) >= 15))) {
  te <- dd[dd$centre == cc, ]; tr <- dd[dd$centre != cc, ]; n <- nrow(te)
  pr <- fit_prior(tr, nut); la <- log(te$soil_loss)
  pub <- dm_nutrient_empirical(te$soil_loss, "N", "soil"); pub <- if (is.list(pub)) pub$mean else pub
  # (a) published equation, plug-in;  (b) hierarchical model, no local data
  ra <- regret(pub > Tn, te$N)
  rb <- regret(decide(pr$b * la + pr$m, pr$tau2 + pr$s2), te$N)
  # (c) k random local plots;  (d) k plots chosen by value of information
  upd <- function(cal) { post <- posterior(pr, log(te$N[cal]) - pr$b * la[cal])
    regret(decide(pr$b * la + post[1], post[2] + pr$s2, te$N[cal], cal), te$N) }
  rc <- mean(replicate(R, upd(sample(n, k))))
  Sig <- matrix(pr$tau2, n, n) + diag(pr$s2, n)
  pred <- dm_gaussian(pr$b * la + pr$m, cov = Sig)
  des <- dm_design(pred, loss_log, n = k, noise = 1e-8)
  rd <- upd(des$site)
  voi_pred <- sum(des$voi)
  res <- rbind(res, data.frame(centre = cc, plots = n, exceed = mean(te$N > Tn),
                               published = ra, model_no_local = rb, random5 = rc, voi5 = rd))
}
res$cost_units <- NULL
print(res, digits = 3)
w <- res$plots / sum(res$plots)
cat("\nPlot-weighted mean regret per plot:\n"); print(round(colSums(res[, 4:7] * w), 3))
saveRDS(res, file.path(outd, "out", "decision.rds"))

# Curve: realised regret (random local plots) versus expected loss predicted by the VOI design
ks <- c(0, 1, 2, 3, 5, 8, 10); curve <- NULL
for (cc in res$centre) {
  te <- dd[dd$centre == cc, ]; tr <- dd[dd$centre != cc, ]; n <- nrow(te)
  pr <- fit_prior(tr, nut); la <- log(te$soil_loss)
  Sig <- matrix(pr$tau2, n, n) + diag(pr$s2, n)
  des <- dm_design(dm_gaussian(pr$b * la + pr$m, cov = Sig), loss_log, n = max(ks), noise = 1e-8)
  expl <- attr(des, "remaining_loss") / n
  for (k in ks) {
    real <- if (k == 0) regret(decide(pr$b * la + pr$m, pr$tau2 + pr$s2), te$N) else
      mean(replicate(R, { cal <- sample(n, k); post <- posterior(pr, log(te$N[cal]) - pr$b * la[cal])
        regret(decide(pr$b * la + post[1], post[2] + pr$s2, te$N[cal], cal), te$N) }))
    curve <- rbind(curve, data.frame(centre = cc, n = n, k = k, realised = real, expected = expl[k + 1]))
  }
}
agg <- do.call(rbind, lapply(split(curve, curve$k), function(s)
  data.frame(k = s$k[1], realised = weighted.mean(s$realised, s$n), expected = weighted.mean(s$expected, s$n))))
print(agg, digits = 3)
saveRDS(list(curve = curve, agg = agg), file.path(outd, "out", "voi_curve.rds"))
x <- readRDS(file.path(outd, "out", "nutrients.rds")); d <- x$d
cal <- readRDS(file.path(outd, "out", "calibration_curve.rds")); vc <- readRDS(file.path(outd, "out", "voi_curve.rds"))$agg
dec <- readRDS(file.path(outd, "out", "decision.rds"))
names_c <- c(BL = "Bellary", CH = "Chandigarh", DN = "Dehradun", DT = "Datia", KR = "Koraput", KT = "Kota", OT = "Ooty")
cols <- c(BL = "#1b9e77", CH = "#999999", DN = "#d95f02", DT = "#7570b3", KR = "#e7298a", KT = "#66a61e", OT = "#1f78b4")
png(file.path(outd, "out", "fig3.png"), 3300, 1050, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1), mgp = c(2.5, 0.8, 0), cex = 0.9)
ok <- d$soil_loss > 0.5 & d$N > 0 & !is.na(d$N)
plot(d$soil_loss[ok], d$N[ok], log = "xy", pch = 16, cex = 0.7, col = adjustcolor(cols[d$centre[ok]], 0.7),
     xlab = expression("Soil loss (t ha"^-1*")"), ylab = expression("N loss (kg ha"^-1*")"), main = "(a) N loss by centre")
for (cc in names(cols)) { s <- d[ok & d$centre == cc, ]; if (nrow(s) < 8) next
  f <- lm(log(N) ~ log(soil_loss), s); xr <- range(s$soil_loss); xs <- exp(seq(log(xr[1]), log(xr[2]), length = 20))
  lines(xs, exp(predict(f, data.frame(soil_loss = xs))), col = cols[cc], lwd = 2) }
present <- names(cols)[names(cols) %in% d$centre[ok]]
legend("bottomright", names_c[present], col = cols[present], pch = 16, bty = "n", cex = 0.75)
plot(NA, xlim = c(0, 10), ylim = c(-0.45, 0.55), xlab = "Local calibration plots", ylab = "NSE (log scale)",
     main = "(b) Prediction in a new centre")
abline(h = 0, col = "grey60", lty = 2)
nc <- c(N = "#d95f02", P = "#7570b3", K = "#1b9e77")
for (nu in names(nc)) { s <- cal[cal$nutrient == nu, ]; lines(s$k, s$NSE_log, type = "b", pch = 16, col = nc[nu], lwd = 2) }
legend("bottomright", c("N", "P", "K"), col = nc, lwd = 2, pch = 16, bty = "n")
plot(vc$k, vc$realised, type = "b", pch = 16, lwd = 2, col = "#d95f02", ylim = c(0, 0.75),
     xlab = "Local calibration plots", ylab = "Regret per plot", main = "(c) Decision: N loss > 10 kg/ha")
lines(vc$k, vc$expected, type = "b", pch = 1, lwd = 2, lty = 2, col = "#d95f02")
w <- dec$plots / sum(dec$plots); abline(h = sum(dec$published * w), col = "grey40", lty = 3, lwd = 2)
legend("topright", c("realised", "expected (TerravAI)", "published equation, in-sample"), col = c("#d95f02", "#d95f02", "grey40"),
       lty = c(1, 2, 3), pch = c(16, 1, NA), lwd = 2, bty = "n", cex = 0.8)
dev.off()
