# Reproduces Section 5.5 of the TerravAI paper: RUSLE on global cropland runoff plots.
# Data: SoilErosionDB (Jian et al.), https://github.com/jinshijian/SoilErosionDB, file data/AWESOME_V1.csv
# Usage: TERRAVAI_SEDB=AWESOME_V1.csv TERRAVAI_OUT=outdir Rscript global_analysis.R
sedb <- Sys.getenv("TERRAVAI_SEDB", "AWESOME_V1.csv"); outd <- Sys.getenv("TERRAVAI_OUT", tempdir()); dir.create(file.path(outd, "out"), showWarnings = FALSE)
raw <- read.csv(text = iconv(readLines(sedb, warn = FALSE), "latin1", "UTF-8"), stringsAsFactors = FALSE)
suppressMessages(library(terravai))
db <- raw
p <- subset(db, Meas_method == "Runoff plot" & IGBP == "CRO" & !is.na(ER_annual) & !is.na(Slope) &
              !is.na(L) & !is.na(MAP) & MAP > 0 & L > 0 & Slope >= 0 & Slope < 60)
p <- p[is.na(p$Duplicate_record) | p$Duplicate_record == "", ]
cat("cropland runoff-plot records:", nrow(p), " sites:", length(unique(p$Site_ID)), " studies:",
    length(unique(p$Study_number)), " countries:", length(unique(p$Country)), "\n")
print(sort(table(p$Country), decreasing = TRUE)[1:12])
# Erosivity from mean annual precipitation (Renard and Freimund, 1994)
Rm <- ifelse(p$MAP < 850, 0.04830 * p$MAP^1.610, 587.8 - 1.219 * p$MAP + 0.004105 * p$MAP^2)
# Erodibility: EPIC where texture and SOC are reported, else a generic prior
has_tex <- !is.na(p$Soil_sand) & !is.na(p$Soil_silt) & !is.na(p$Soil_clay) & !is.na(p$Soil_SOC)
Km <- rep(0.03, nrow(p)); Kcv <- rep(0.5, nrow(p))
Km[has_tex] <- dm_k_epic(pmin(pmax(p$Soil_sand[has_tex], 0), 100), pmin(pmax(p$Soil_silt[has_tex], 0), 100),
                         pmin(pmax(p$Soil_clay[has_tex], 0), 100), pmax(p$Soil_SOC[has_tex], 0))
Kcv[has_tex] <- 0.3
cat("records with texture-based K:", sum(has_tex), "\n")
LSs <- suppressWarnings(dm_ls_factor(p$Slope, length = p$L, method = "wischmeier"))
LSn <- suppressWarnings(dm_ls_factor(p$Slope, length = p$L, method = "wischmeier", steepness = "nearing"))
obs <- p$ER_annual
run <- function(LS, seed) dm_rusle(R = dm_uncertain(Rm, cv = 0.4, dist = "lognormal"),
                                   K = dm_uncertain(Km, cv = Kcv, dist = "lognormal"),
                                   LS = LS, C = dm_uncertain(rep(0.233, nrow(p)), cv = 0.6, dist = "lognormal"),
                                   P = 1, n_draws = 2000, seed = seed)
A_std <- run(LSs, 1); A_near <- run(LSn, 1)
cat("\n--- RUSLE, standard slope factor\n"); print(dm_validate(obs, A_std))
cat("\n--- RUSLE, Nearing slope factor\n"); print(dm_validate(obs, A_near))
cat("\nlog-scale NSE std / nearing:",
    round(1 - sum((log(A_std$mean) - log(pmax(obs, 1e-3)))^2) / sum((log(pmax(obs,1e-3)) - mean(log(pmax(obs,1e-3))))^2), 3), "/",
    round(1 - sum((log(A_near$mean) - log(pmax(obs, 1e-3)))^2) / sum((log(pmax(obs,1e-3)) - mean(log(pmax(obs,1e-3))))^2), 3), "\n")
cat("median obs / median pred:", round(median(obs), 2), round(median(A_std$mean), 2), "\n")
cat("Spearman obs vs pred:", round(cor(obs, A_std$mean, method = "spearman"), 3), "\n")
print(dm_convergence(A_std, dm_loss_tolerance(T_limit = 11.2)))
for (Tl in c(5, 11.2)) {
  cat("\n=== decision: soil loss >", Tl, "t/ha/yr; share exceeding =", round(mean(obs > Tl), 3), "\n")
  ev <- dm_evaluate(obs, list(decision_rule = A_std, plug_in = A_std$mean), dm_loss_tolerance(T_limit = Tl), ratios = c(1, 3, 10))
  print(as.data.frame(unclass(ev))[, c("method", "ratio", "regret", "missed", "unneeded")], digits = 3)
}
ind <- p$Country == "India"
cat("\nIndia records:", sum(ind), "\n"); if (sum(ind) >= 5) print(dm_validate(obs[ind], dm_draws(A_std$draws[ind, , drop = FALSE])))
saveRDS(list(p = p, Rm = Rm, Km = Km, A_std = A_std, A_near = A_near), file.path(outd, "out", "global.rds"))
suppressMessages(library(terravai))
g <- readRDS(file.path(outd, "out", "global.rds")); p <- g$p; A <- g$A_std
p$pred <- A$mean; p$obs <- p$ER_annual
# slope-unit sensitivity: if Slope were percent rather than degrees
deg_from_pct <- atan(p$Slope / 100) * 180 / pi
LSpct <- suppressWarnings(dm_ls_factor(deg_from_pct, length = p$L, method = "wischmeier"))
cat("median pred if slope in %:", round(median(g$Rm * g$Km * LSpct * 0.233), 2), " (degrees:", round(median(p$pred), 2), "; obs", round(median(p$obs), 2), ")\n")
cat("zero observations:", sum(p$obs <= 0), "\n")
studies <- names(which(table(p$Study_number) >= 15)); cat("studies with >= 15 records:", length(studies), " records:", sum(p$Study_number %in% studies), "\n")
Tl <- 11.2; r <- 3
regret <- function(act, y) mean(ifelse(act & y <= Tl, 1, 0) + ifelse(!act & y > Tl, r, 0))
lnse <- function(o, pm) { o <- pmax(o, 0.01); 1 - sum((log(pm) - log(o))^2) / sum((log(o) - mean(log(o)))^2) }
set.seed(5); ks <- c(0, 3, 5, 10); R <- 30
res <- NULL
for (st in studies) {
  te <- p[p$Study_number == st, ]; tr <- p[p$Study_number != st, ]; n <- nrow(te)
  D_raw <- A$draws[match(rownames(te), rownames(p)), , drop = FALSE]
  base <- data.frame(study = st, n = n, k = NA, method = c("RUSLE plug-in", "RUSLE decision rule"),
                     lnse = lnse(te$obs, te$pred),
                     cov90 = mean(te$obs >= apply(D_raw, 1, quantile, .05) & te$obs <= apply(D_raw, 1, quantile, .95)),
                     regret = c(regret(te$pred > Tl, te$obs), regret(rowMeans(D_raw > Tl) > 1 / (1 + r), te$obs)))
  res <- rbind(res, base)
  for (k in ks) for (rep in seq_len(if (k == 0) 1 else R)) {
    cal <- if (k == 0) integer(0) else sample(which(te$obs > 0), min(k, sum(te$obs > 0)))
    f <- dm_transfer(tr$pred, tr$obs, tr$Study_number, new_x = te$pred,
                     local_x = if (k) te$pred[cal] else NULL, local_y = if (k) te$obs[cal] else NULL,
                     n_draws = 1000, seed = rep)
    D <- f$draws; med <- exp(attr(f, "log_mean"))
    act <- rowMeans(D > Tl) > 1 / (1 + r); if (k) act[cal] <- te$obs[cal] > Tl
    ev <- setdiff(seq_len(n), cal)
    res <- rbind(res, data.frame(study = st, n = n, k = k, method = "calibrated",
                                 lnse = lnse(te$obs[ev], med[ev]),
                                 cov90 = mean(te$obs[ev] >= apply(D[ev, , drop = FALSE], 1, quantile, .05) &
                                                te$obs[ev] <= apply(D[ev, , drop = FALSE], 1, quantile, .95)),
                                 regret = regret(act, te$obs)))
  }
}
agg <- aggregate(cbind(lnse, cov90, regret) ~ method + k, data = transform(res, k = ifelse(is.na(k), -1, k)),
                 FUN = mean)
# weight studies equally after averaging repetitions
agg2 <- do.call(rbind, lapply(split(res, list(res$method, ifelse(is.na(res$k), -1, res$k)), drop = TRUE), function(s) {
  ps <- aggregate(cbind(lnse, cov90, regret) ~ study + n, s, mean)
  data.frame(method = s$method[1], k = s$k[1], studies = nrow(ps),
             lnse_median = median(ps$lnse), cov90 = weighted.mean(ps$cov90, ps$n), regret = weighted.mean(ps$regret, ps$n)) }))
print(agg2[order(agg2$method, agg2$k), ], digits = 3, row.names = FALSE)
f_all <- dm_transfer(p$pred, p$obs, p$Study_number, new_x = 10)
print(round(attr(f_all, "params"), 3))
saveRDS(list(res = res, agg = agg2, params = attr(f_all, "params")), file.path(outd, "out", "global_calib.rds"))
suppressMessages(library(terravai))
g <- readRDS(file.path(outd, "out", "global.rds")); p <- g$p; A <- g$A_std
cal <- readRDS(file.path(outd, "out", "global_calib.rds")); a <- cal$agg
fa <- dm_transfer(A$mean, p$ER_annual, p$Study_number, new_x = A$mean)   # pooled calibration (in-sample, display only)
png(file.path(outd, "fig4.png"), 3300, 1050, res = 300)
par(mfrow = c(1, 3), mar = c(4.2, 4.4, 2.6, 1), mgp = c(2.5, 0.8, 0), cex = 0.9)
o <- pmax(p$ER_annual, 0.01); tk <- c(0.01, 0.1, 1, 10, 100, 1000)
plot(A$mean, o, log = "xy", pch = 16, cex = 0.5, col = adjustcolor("grey30", 0.35), xlim = c(0.01, 3000), ylim = c(0.01, 3000),
     axes = FALSE, xlab = expression("RUSLE prediction (t ha"^-1*" yr"^-1*")"), ylab = expression("Measured soil loss (t ha"^-1*" yr"^-1*")"),
     main = "(a) RUSLE on 864 cropland plots")
axis(1, tk, format(tk, scientific = FALSE, drop0trailing = TRUE)); axis(2, tk, format(tk, scientific = FALSE, drop0trailing = TRUE)); box()
abline(0, 1, lty = 2, col = "grey50"); xs <- exp(seq(log(0.05), log(3000), length = 50))
b <- cal$params[["b"]]; m <- cal$params[["m"]]
lines(xs, exp(m + b * log(xs)), col = "#1f78b4", lwd = 2.5)
legend("topleft", c("1:1", "calibrated across studies"), lty = c(2, 1), col = c("grey50", "#1f78b4"), lwd = c(1, 2.5), bty = "n", cex = 0.8)
cc <- a[a$method == "calibrated", ]; cc <- cc[order(cc$k), ]
raw <- a[a$method == "RUSLE decision rule", ]
plot(cc$k, cc$regret, type = "b", pch = 16, lwd = 2, col = "#d95f02", ylim = c(0, 0.75),
     xlab = "Local calibration plots", ylab = "Regret per plot", main = "(b) Decision: soil loss > 11.2 t/ha")
abline(h = raw$regret, lty = 3, lwd = 2, col = "grey40")
legend("bottomleft", c("calibrated", "uncalibrated RUSLE"), col = c("#d95f02", "grey40"), lty = c(1, 3), pch = c(16, NA), lwd = 2, bty = "n", cex = 0.8)
plot(cc$k, 100 * cc$cov90, type = "b", pch = 16, lwd = 2, col = "#1f78b4", ylim = c(0, 100),
     xlab = "Local calibration plots", ylab = "Coverage of 90% interval (%)", main = "(c) Honesty of uncertainty")
abline(h = 90, lty = 2, col = "grey50"); abline(h = 100 * raw$cov90, lty = 3, lwd = 2, col = "grey40")
legend("bottomleft", c("calibrated", "uncalibrated RUSLE", "nominal 90%"), col = c("#1f78b4", "grey40", "grey50"), lty = c(1, 3, 2), pch = c(16, NA, NA), lwd = c(2, 2, 1), bty = "n", cex = 0.8)
dev.off()
