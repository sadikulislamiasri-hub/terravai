# terravai 0.1.0

* First release of TerravAI (package `terravai`).
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
