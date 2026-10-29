# terravai 0.1.0

Authors: Sadikul Islam (maintainer), Rajesh Kaushal, Shakir Ali, M Shankar,
M Muruganandam and Rama Pal.

* First release of TerravAI (Threshold Expected-Regret and Risk Analysis with
  Value of information and Artificial Intelligence), package `terravai`.
* `dm_transfer()`: calibrates any formula (for example RUSLE) against
  measurements across sites, with local calibration at a new site.
* `dm_nutrient_transfer()`: multi-site nutrient-loss model (shared slope,
  site effect) calibrated with a few local plots; predictive draws share the
  site effect, so value of information can price local calibration.
* Second review round:
  * `dm_fill_sinks()` rewritten in C (priority-flood, Barnes et al. 2014):
    same surface as before, about 100 times faster (4 million cells in about
    1 s).
  * `steepness = "nearing"` in `dm_ls_factor()` and `dm_erosion_raster()`:
    Nearing (1997) slope factor for steep terrain.
  * `dm_gully_risk()` and `SPI`/`TWI` layers in `dm_erosion_raster()` to
    screen gully and ravine risk.
  * `dm_krige(nmax = )`: local kriging for maps of any size.
  * `dm_validate()`: NSE, KGE, PBIAS, RMSE, R2, CRPS and interval coverage
    against runoff-plot measurements.
  * `dm_erosivity_table()` and `dm_erosivity_india(equation = )`: regional
    erosivity relations (13 Bundelkhand districts, Gupta et al. 2024), with a
    warning outside the fitted rainfall range.
* Correctness review (before first release):
  * `dm_fill_sinks()`: Planchon-Darboux pit filling; `dm_erosion_raster()`
    now fills sinks by default (`fill_sinks = TRUE`). On a noisy test DEM this
    cut false outlets from 4,217 to 162.
  * `dm_runoff_cn_daily()`: daily curve number runoff with antecedent
    moisture from the previous five days' rainfall.
  * `dm_flow_from_dem()` uses true cell spacing (`dx`, `dy`), including
    metres on longitude/latitude grids.
  * Formulas now keep the spatial correlation of Gaussian predictive inputs
    with a full covariance (previously sampled independently).
  * Formula inputs of incompatible lengths are now an error instead of being
    silently recycled.
  * `dm_outlet_load()` routes Monte Carlo draws one at a time in linear time
    instead of building a dense routing matrix.
  * Warnings for steep slopes beyond the (R)USLE plot data, and for ensemble
    value of information with fewer than 1,000 draws.
* Runoff, soil loss and nutrient loss:
  * `dm_runoff_cn()`, `dm_cn_amc()`, `dm_peak_rational()`: curve-number runoff,
    antecedent moisture adjustment and rational-method peak flow.
  * `dm_rusle()` with `dm_ls_factor()`, `dm_c_ndvi()`, `dm_k_epic()` and
    `dm_fournier()`; `dm_musle()`; `dm_sediment_yield()`.
  * `dm_nutrient_sediment()`, `dm_nutrient_runoff()`, `dm_enrichment()`:
    nutrient loss by mass balance.
  * `dm_nutrient_empirical()` and `dm_nutrient_coef()`: empirical N, P and K
    loss equations from ICAR-IISWC runoff-plot data (first released in
    SoilRunNutri; contributors Soumya Pant, Gopal Kumar and M. Madhu).
  * `dm_uncertain()`: uncertain inputs carried through every formula.
  * `dm_flow_accumulation()` (linear time) and `dm_outlet_load()`.
  * `dm_loss_tolerance()` and `dm_loss_nutrient()` decision presets.
  * `dm_simulate_watershed()`: example micro-watershed with known truth.
* Indian standards:
  * `dm_psd_isss()`, `dm_psd_is1498()`, `dm_psd_to_usda()`: convert
    International (ISSS), Indian Standard IS 1498 or any particle-size data to
    USDA fractions by log-linear interpolation.
  * `dm_shc_table()`, `dm_shc_rating()`: Soil Health Card ratings (ICAR-IISS
    limits), with class probabilities for uncertain maps.
  * `dm_erosivity_india()`: erosivity from annual rainfall (Singh, Chandra and
    Babu, 1981).
  * `dm_slope_class_india()`: AIS&LUS slope classes A-I.
  * `dm_erosion_class()`: erosion severity classes, default E0-E4 or custom.
* Uncertainty and convergence:
  * `dm_uncertain(correlation = )`: input errors can be independent, common
    to all locations (e.g. a regional coefficient) or spatially correlated
    (Gaussian copula); per-location distributions are unchanged.
  * `dm_convergence()`: Monte Carlo standard errors, running estimates,
    decisions still within Monte Carlo noise, and the number of draws needed.
  * Closed-form value of information now integrates separately on each side
    of the decision switch point: 50-100 times more accurate than before at
    the same cost (below 0.1% error with the default 20 nodes per side) and
    checked against Monte Carlo.
  * `dm_fit_cov()`: multiple starting values, and warnings when the optimiser
    does not converge or a parameter ends at its search limit.
  * `dm_krige()`: handles a single prediction location and near-singular
    covariance matrices (with a warning).
* Use anywhere (pan-India and global):
  * `dm_erosion_raster()`: RUSLE chain on rasters for any region, on
    projected or longitude/latitude grids, with SoilGrids units supported.
  * `dm_cn_table()` and `dm_cn_lookup()` (USDA TR-55 curve numbers),
    `dm_texture_class()` (USDA texture classes), `dm_hsg_texture()`,
    `dm_p_factor()` (Wischmeier and Smith support-practice factors).
  * Vectorised D8 flow directions and linear-time flow accumulation (about a
    million cells per second).
  * New vignette on runoff, soil loss and nutrient loss.
* Predictive distributions from means and standard deviations, covariance
  matrices, simulation draws or `terra` rasters.
* AI and machine learning:
  * `dm_fit_qrf()`: quantile regression forests ('ranger') with a spatial
    Gaussian copula estimated from out-of-bag residuals.
  * `dm_quantiles()`: predictive quantiles from any model, joined across space
    with an optional spatial copula.
  * `dm_conformal()`: split-conformal predictive distributions for any point
    predictor, with optional local scaling.
  * `dm_voi(method = "ensemble")`: value of information for non-Gaussian
    ensembles by an ensemble Kalman preposterior update with localisation.
  * `dm_narrate()`: plain-language summaries of decision maps.
* Threshold decisions with fixed or linear losses defined by a single cost ratio.
* Decision maps: optimal action, probability of a wrong decision, expected loss
  and break-even cost ratio.
* Closed-form value-of-information maps under Gaussian updating.
* Greedy budget-constrained sampling designs with lazy evaluation and a
  value-based stopping rule.
* Decision-based evaluation of competing maps, ranking reversals and decision
  curves across thresholds.
* Scale operators for districts (block means) and watersheds (D8 flow
  accumulation).
* A dependency-free kriging engine with (RE)ML covariance fitting.
