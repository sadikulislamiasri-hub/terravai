# TerravAI <img src="man/figures/logo.png" align="right" height="139" alt="TerravAI logo" />

<!-- badges: start -->
[![R-CMD-check](https://github.com/GITHUB_USER/terravai/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/GITHUB_USER/terravai/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

**TerravAI: Threshold Expected-Regret and Risk Analysis with Value of information and Artificial Intelligence**

*AI-driven runoff, soil and nutrient loss decisions*

TerravAI (R package `terravai`) quantifies **runoff, soil loss and nutrient
loss** with their uncertainty, from formulas (curve number, RUSLE, MUSLE,
nutrient mass balance, ICAR-IISWC N/P/K equations) or AI models, and turns
them into conservation decisions. It turns the predictive distribution of any
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

## Runoff, soil loss and nutrient loss

```r
library(terravai)
ws <- dm_simulate_watershed(seed = 3)          # example micro-watershed
g  <- ws$grid
xy <- g[, c("x", "y")]

# Soil loss by RUSLE with uncertain erosivity and cover
soil <- dm_rusle(R  = dm_uncertain(g$R, cv = 0.15, lower = 0, coords = xy),
                 K  = g$K, LS = g$ls,
                 C  = dm_uncertain(g$C, sd = 0.05, lower = 0, upper = 1, coords = xy),
                 seed = 1)

# Where is soil loss above the tolerable limit (11.2 t/ha/yr)?
dec <- dm_decide(soil, dm_loss_tolerance())
plot(dec, "act"); dm_narrate(dec, area = "the micro-watershed")

# Nutrient loss: mass balance and ICAR-IISWC empirical equations
N_mb  <- dm_nutrient_sediment(soil, content = 0.08, unit = "percent", enrichment = 1.6)
N_emp <- dm_nutrient_empirical(soil, "N", "soil")

# Load at the outlet and the best places for new runoff plots
dm_outlet_load(soil$mean, ws$to, delivery = 0.9, cell_area = 0.09)
plot(dm_voi(soil, dm_loss_tolerance(), method = "ensemble", seed = 1))
```

## Indian standards

| Need | Function | Basis |
|---|---|---|
| Convert laboratory particle-size data to USDA fractions | `dm_psd_isss()`, `dm_psd_is1498()`, `dm_psd_to_usda()` | International (ISSS) system, Indian Standard IS 1498, any cumulative distribution |
| USDA texture class and hydrologic soil group | `dm_texture_class()`, `dm_hsg_texture()` | USDA texture triangle |
| Soil Health Card ratings (OC, N, P, K, S, Zn, Cu, Fe, Mn, B), with class probabilities for uncertain maps | `dm_shc_rating()`, `dm_shc_table()` | ICAR-IISS rating limits |
| Rainfall erosivity from annual rainfall | `dm_erosivity_india()` | R = 79 + 0.363 P (Singh, Chandra and Babu, 1981) |
| Slope classes A-I | `dm_slope_class_india()` | All India Soil and Land Use Survey |
| Erosion severity classes E0-E4 (or your own scheme) | `dm_erosion_class()` | USLE-based erosion mapping in India |
| Empirical N, P, K loss from runoff plots | `dm_nutrient_empirical()` | ICAR-IISWC experiments |
| Site-specific tolerance limits | `dm_loss_tolerance()` | e.g. 2.5-12.5 t/ha/yr estimated for Indian mapping units |

```r
usda <- dm_psd_isss(clay = 28, silt = 14, fine_sand = 38, coarse_sand = 20)
dm_texture_class(usda$sand, usda$silt, usda$clay)        # "clay loam"
dm_shc_rating(c(180, 350, 620), "N")                       # low, medium, high
R <- dm_erosivity_india(1400)                              # from annual rainfall
dm_erosion_class(dm_rusle(R, 0.03, 3, 0.25))
```

## Uncertainty you can trust

* **Input uncertainty:** any input can be uncertain with `dm_uncertain()`, as
  independent errors, an error common to the whole area (such as a regional
  equation's coefficients) or spatially correlated errors.
* **Convergence checks:** `dm_convergence()` reports Monte Carlo standard
  errors, flags decisions that could still change by chance, and recommends
  the number of draws.
* **Verified numerics:**
  * Monte Carlo errors fall at the theoretical 1/sqrt(n) rate.
  * The closed-form value of information agrees with Monte Carlo and is
    accurate to better than 0.1%.
  * The ensemble method converges to the exact answer as draws increase.
  * Covariance fitting uses several starting values and warns when a
    parameter is poorly determined.

## Use anywhere: pan-India and global

Every tool runs on any area with freely available data, on projected or
longitude/latitude grids:

```r
library(terra)
dem  <- rast("dem.tif")                 # e.g. SRTM / Copernicus DEM, or Bhuvan CartoDEM
ndvi <- rast("ndvi.tif")                # e.g. Sentinel-2 or MODIS NDVI
sand <- rast("sand.tif"); silt <- rast("silt.tif"); clay <- rast("clay.tif")
soc  <- rast("soc.tif")                 # e.g. SoilGrids (g/kg, dg/kg) or national soil maps
R    <- rast("erosivity.tif")           # national or global erosivity map

out <- dm_erosion_raster(dem, R, ndvi = ndvi, sand = sand, silt = silt, clay = clay,
                         soc = soc, texture_unit = "g_kg", soc_unit = "dg_kg")
plot(out$A)                             # soil loss, t/ha/yr
writeRaster(out, "soil_loss.tif")
```

| Need | Tools |
|---|---|
| Curve numbers for any land use and soil | `dm_cn_table()` (USDA TR-55), `dm_cn_lookup()` |
| Hydrologic soil group where only texture is known | `dm_texture_class()` (USDA), `dm_hsg_texture()` |
| Conservation practice factor | `dm_p_factor()` (contouring, strip cropping) |
| Region-specific limits | `dm_loss_tolerance(T_limit = )`, or a site-specific tolerance map (see `?dm_loss_tolerance`) |
| Large rasters | D8 routing and flow accumulation handle about a million cells in about a second |

Regional coefficients that vary between places (rainfall erosivity relations,
enrichment ratios, local C, P and curve-number tables, tolerance limits) are
user inputs rather than hidden defaults, so the package applies equally to
any state in India or any country.

## Quick start (general decision engine)

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
| **Runoff:** curve number, moisture adjustment, rational peak flow | `dm_runoff_cn()`, `dm_cn_amc()`, `dm_peak_rational()` |
| **Soil loss:** RUSLE and factors, MUSLE, sediment delivery | `dm_rusle()`, `dm_ls_factor()`, `dm_c_ndvi()`, `dm_k_epic()`, `dm_fournier()`, `dm_musle()`, `dm_sediment_yield()` |
| **Nutrient loss:** mass balance, enrichment, ICAR-IISWC N/P/K equations | `dm_nutrient_sediment()`, `dm_nutrient_runoff()`, `dm_enrichment()`, `dm_nutrient_empirical()`, `dm_nutrient_coef()` |
| **Uncertain inputs** for any formula | `dm_uncertain()` |
| **Watershed routing** | `dm_flow_from_dem()`, `dm_flow_accumulation()`, `dm_outlet_load()` |
| **Conservation decisions** | `dm_loss_tolerance()`, `dm_loss_nutrient()` |
| **Example micro-watershed** with known truth | `dm_simulate_watershed()` |
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

## Authors

- **Sadikul Islam** (author, maintainer) [![ORCID](https://img.shields.io/badge/ORCID-0000--0003--2924--7122-a6ce39)](https://orcid.org/0000-0003-2924-7122)
- Rajesh Kaushal (author)
- Shakir Ali (author)
- M Shankar (author)
- M Muruganandam (author)
- Rama Pal (author)

## Contributors

The empirical nitrogen, phosphorus and potassium loss equations
(`dm_nutrient_empirical()`) were first released in the Python package
SoilRunNutri by Soumya Pant, Sadikul Islam, Gopal Kumar and M. Madhu; see
`inst/COPYRIGHTS`.

## Acknowledgements

The authors gratefully acknowledge the Indian Council of Agricultural Research (ICAR), New Delhi, for providing the necessary research support. This work was carried out as part of the institute-funded project “Development of AI-driven solutions for quantifying runoff induced soil loss and nutrient loss” (NRMAIISWCSIL202400300220) at ICAR–Indian Institute of Soil and Water Conservation, Dehradun.

## License

MIT © terravai authors (Sadikul Islam, Rajesh Kaushal, Shakir Ali, M Shankar,
M Muruganandam and Rama Pal)
