# Reproduces every number and figure in the TerravAI software paper.
# Run with: Rscript system.file("paper", "reproduce.R", package = "terravai")
library(terravai)
out_dir <- if (nzchar(Sys.getenv("TERRAVAI_PAPER_DIR"))) Sys.getenv("TERRAVAI_PAPER_DIR") else tempdir()
res <- list()
ws <- dm_simulate_watershed(n_side = 30, seed = 3)
g <- ws$grid; xy <- g[, c("x", "y")]
res$cells <- nrow(g); res$area_ha <- nrow(g) * ws$cell_size^2 / 1e4

# 1. Soil loss with uncertain erosivity and cover
soil <- dm_rusle(R = dm_uncertain(g$R, cv = 0.15, lower = 0, coords = xy), K = g$K, LS = g$ls,
                 C = dm_uncertain(g$C, sd = 0.05, lower = 0, upper = 1, coords = xy),
                 n_draws = 1000, seed = 2)
res$soil_mean <- mean(soil$mean); res$soil_range <- range(soil$mean)
res$true_mean <- mean(g$soil_loss)

# 2. Convergence
cv <- dm_convergence(soil, dm_loss_tolerance())
res$conv <- cv

# 3. Decision
loss <- dm_loss_tolerance(T_limit = 11.2, ratio = 3)
dec <- dm_decide(soil, loss)
res$share_act <- mean(dec$act); res$share_uncertain <- mean(dec$p_wrong > 0.2)
res$true_share <- mean(g$soil_loss > 11.2)

# 4. Accuracy versus decision quality
ev <- dm_evaluate(g$soil_loss, list(decision_rule = soil, plug_in = soil$mean), dm_loss_tolerance(),
                  ratios = c(1, 3, 10))
res$eval <- summary(ev)

# 5. Validation
res$valid <- dm_validate(g$soil_loss, soil)

# 6. Correlated erosivity error and watershed total
tot_sd <- function(corr) {
  s <- dm_rusle(R = dm_uncertain(g$R, cv = 0.15, lower = 0, coords = xy, correlation = corr),
                K = g$K, LS = g$ls, C = g$C, n_draws = 1000, seed = 5)
  t <- colSums(s$draws); c(mean = mean(t), sd = sd(t))
}
res$total_indep <- tot_sd("independent"); res$total_common <- tot_sd("common")

# 7. Value of information and plot design
cand <- seq(1, nrow(g), by = 3)
v <- suppressWarnings(dm_voi(soil, dm_loss_tolerance(), method = "ensemble", candidates = cand, seed = 4))
res$voi_total <- sum(v$voi, na.rm = TRUE); res$voi_top <- head(sort(v$voi, decreasing = TRUE), 5)
res$voi_obj <- v

# 8. Closed-form VOI accuracy against Monte Carlo
mu <- 10; s <- 3; tau <- 2; Tl <- 11.2; r <- 3
set.seed(1); M <- 2e6
th <- rnorm(M, mu, s); yobs <- th + rnorm(M, 0, tau)
w <- s^2 / (s^2 + tau^2); mpost <- mu + w * (yobs - mu); spost <- sqrt(s^2 * (1 - w))
pexc <- 1 - pnorm((Tl - mpost) / spost)
prior_loss <- min(r * (1 - pnorm((Tl - mu) / s)), pnorm((Tl - mu) / s))
post_loss <- mean(pmin(r * pexc, 1 - pexc))
res$voi_mc <- prior_loss - post_loss

# 9. Daily versus annual-shortcut curve-number runoff
set.seed(7)
rain <- ifelse(runif(122) < 0.35, rexp(122, 1 / 28), 0)
daily <- dm_runoff_cn_daily(rain, CN = c(70, 80, 88))
shortcut <- sapply(c(70, 80, 88), function(cn) {
  S <- 25400 / cn - 254; P <- sum(rain) / 25; 25 * ifelse(P > 0.2 * S, (P - 0.2 * S)^2 / (P + 0.8 * S), 0)
})
res$rain_total <- sum(rain); res$daily <- daily$total; res$shortcut <- shortcut

# 10. Scaling benchmarks
bench <- t(sapply(c(500, 1000, 2000), function(n) {
  set.seed(1)
  z <- outer(1:n, 1:n, function(i, j) 0.5 * i + 0.3 * j) + matrix(rnorm(n * n, 0, 2), n)
  t1 <- system.time(f <- dm_fill_sinks(z))[[3]]
  t2 <- system.time(to <- dm_flow_from_dem(f))[[3]]
  t3 <- system.time(acc <- dm_flow_accumulation(to))[[3]]
  pits <- sum(is.na(dm_flow_from_dem(z))); outs <- sum(is.na(to))
  c(cells = n * n, fill_s = t1, d8_s = t2, acc_s = t3, sinks_before = pits, sinks_after = outs)
}))
res$bench <- bench
sim <- dm_simulate(n_side = 40, n_obs = 300, seed = 5)
o <- sim$obs; gx <- sim$grid[, c("x", "y")]
t_glob <- system.time(f1 <- dm_krige(o[, c("x", "y")], o$z, gx, sim$model, full_cov = FALSE))[[3]]
t_loc <- system.time(f2 <- suppressMessages(dm_krige(o[, c("x", "y")], o$z, gx, sim$model, nmax = 40)))[[3]]
res$krige <- c(cells = nrow(gx), cor_mean = cor(f1$mean, f2$mean), max_sd_diff = max(abs(f1$sd - f2$sd)),
               ms_per_cell = 1000 * t_loc / nrow(gx))

# Figure: decision maps
png(file.path(out_dir, "fig_decision.png"), width = 1800, height = 620, res = 200)
par(mfrow = c(1, 3), mar = c(1, 1, 3, 4))
plot(dec, "act"); plot(dec, "p_wrong"); plot(v, main = "Value of a new runoff plot")
dev.off()
saveRDS(res, file.path(out_dir, "results.rds"))
str(res[setdiff(names(res), c("voi_obj", "conv"))])
print(res$conv); print(res$eval); print(res$valid)
