# terravai 0.1.0

Authors: Sadikul Islam (maintainer), Rajesh Kaushal, Shakir Ali, M Shankar,
M Muruganandam and Rama Pal.

* First release of TerravAI (Threshold Expected-Regret and Risk Analysis with
  Value of information and Artificial Intelligence), package `terravai`.
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
