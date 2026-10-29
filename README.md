# TerravAI <img src="man/figures/logo.png" align="right" height="139" alt="TerravAI logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/GITHUB_USER/terravai/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GITHUB_USER/terravai/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**Value of Information for Spatial Decisions from AI Predictions**

TerravAI (R package `terravai`) turns the predictive distribution of any
spatial model, from kriging to AI and machine-learning ensembles, into
decision maps, expected-loss and regret maps, value-of-information maps and
budget-optimal sampling plans at pixel, district and watershed scales.

Spatial maps are usually chosen and judged by accuracy (RMSE, R²). But maps
are used to *decide*: lime a field, target a micronutrient intervention,
declare a well unsafe. For such threshold decisions, only errors that cross
the threshold matter, and missing a needed action rarely costs the same as an
unneeded one. `terravai` evaluates and designs spatial information for the
decision it informs.

## Installation

```r
# Development version from GitHub
# install.packages("remotes")
remotes::install_github("GITHUB_USER/terravai")
```

## Quick start

```r
library(terravai)

# 1. Predictions from any model (here: built-in kriging on simulated data)
sim <- dm_simulate(n_side = 25, n_obs = 50, seed = 11)
cm  <- dm_fit_cov(sim$obs[, c("x", "y")], sim$obs$z)
fit <- dm_krige(sim$obs[, c("x", "y")], sim$obs$z, sim$grid[, c("x", "y")], model = cm)

# 2. The decision: act above 0.5; a missed action costs 4x an unneeded one
loss <- dm_loss(0.5, "above", ratio = 4, action = "intervene")

# 3. Decision maps
dec <- dm_decide(fit, loss)
plot(dec, "act"); plot(dec, "p_wrong"); plot(dec, "breakeven")

# 4. Where would the next sample be worth most? Plan a survey.
plot(dm_voi(fit, loss, noise = cm$nugget))
plot(dm_design(fit, loss, n = 12, noise = cm$nugget))

# 5. AI predictions: quantile regression forest + spatial copula
sim$grid$ndvi <- sim$grid$z + rnorm(nrow(sim$grid), 0, 0.6)
train <- cbind(sim$obs, ndvi = sim$grid$ndvi[sim$obs$cell])
qrf <- dm_fit_qrf(z ~ ndvi, train, sim$grid, seed = 1)
plot(dm_voi(qrf, loss, noise = cm$nugget, seed = 1))   # ensemble Kalman VOI
dm_narrate(dm_decide(qrf, loss), area = "the district")

# 6. Is the most accurate map the best decision map?
ev <- dm_evaluate(sim$grid$z, list(decision = fit, plugin = fit$mean), loss,
                  ratios = c(1, 2, 4, 8))
summary(ev)
```

## What it does

| Step | Functions |
|---|---|
| Predictive distributions from any model | `dm_gaussian()`, `dm_draws()` (also from `terra` rasters) |
| **AI:** quantile regression forests with a spatial copula | `dm_fit_qrf()`, `dm_quantiles(copula = )` |
| **AI:** conformal predictive distributions for any point predictor | `dm_conformal()` |
| **AI:** value of information for non-Gaussian ensembles (ensemble Kalman) | `dm_voi(method = "ensemble")` |
| **AI:** plain-language summaries of decision maps | `dm_narrate()` |
| Built-in kriging with full predictive covariance | `dm_covmodel()`, `dm_fit_cov()`, `dm_krige()`, `dm_cov()` |
| Threshold decisions with one stakeholder input (cost ratio) | `dm_loss()`, `dm_cutoff()` |
| Decision maps: action, P(wrong), expected loss, break-even ratio | `dm_decide()` |
| Value of information (closed form under Gaussian updating) | `dm_voi()` |
| Budget-constrained sampling with a value-based stopping rule | `dm_design()` |
| Accuracy vs decision regret, ranking reversals, decision curves | `dm_evaluate()`, `dm_threshold_curve()` |
| Districts and watersheds (scale operator `Q = H Z`) | `dm_block_matrix()`, `dm_flow_from_dem()`, `dm_flow_matrix()`, `dm_aggregate()` |
| GIS output | `dm_as_raster()` |
| Simulation with known truth | `dm_simulate()` |

The package has no hard dependencies beyond base R; 'ranger' (forests) and
'terra' (rasters) are optional.

The ensemble value-of-information method is checked against the exact Gaussian
solution: with 1,000 draws the two agree to within about 1% in total and have
a correlation of 0.99 across candidate sites.

## Citation

```r
citation("terravai")
```

## License

MIT © Sadikul Islam
