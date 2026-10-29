#' terravai: AI-Driven Runoff, Soil and Nutrient Loss Decisions
#'
#' TerravAI (Threshold Expected-Regret and Risk Analysis with Value of
#' information and Artificial Intelligence; package `terravai`) turns the predictive distribution of any
#' spatial model, from kriging to AI and machine-learning ensembles, into the
#' quantities a decision-maker needs: which action to take where, the expected
#' loss and probability of a wrong decision, the cost ratio at which each
#' decision would flip, and where new data would be worth the most.
#'
#' @section Workflow:
#' 0. Quantify runoff, soil loss and nutrient loss with [dm_runoff_cn()],
#'    [dm_rusle()] (with [dm_ls_factor()], [dm_c_ndvi()], [dm_k_epic()],
#'    [dm_fournier()]), [dm_musle()], [dm_nutrient_sediment()],
#'    [dm_nutrient_runoff()] and [dm_nutrient_empirical()]; make any input
#'    uncertain with [dm_uncertain()]; route loads with [dm_outlet_load()].
#'    For any region, run the raster chain with [dm_erosion_raster()] and use
#'    the look-up tables [dm_cn_table()], [dm_texture_class()],
#'    [dm_hsg_texture()] and [dm_p_factor()]. Indian standards:
#'    [dm_psd_isss()], [dm_psd_is1498()], [dm_shc_rating()],
#'    [dm_erosivity_india()], [dm_slope_class_india()] and [dm_erosion_class()].
#' 1. Describe predictions with [dm_gaussian()], [dm_draws()],
#'    [dm_quantiles()] or [dm_conformal()], fit a quantile regression forest
#'    with [dm_fit_qrf()], or fit kriging with [dm_fit_cov()] and [dm_krige()].
#' 2. Describe the decision with [dm_loss()], [dm_loss_tolerance()] or
#'    [dm_loss_nutrient()].
#' 3. Map decisions with [dm_decide()] and summarise them in plain language
#'    with [dm_narrate()].
#' 4. Map the value of new samples with [dm_voi()] and plan a survey with
#'    [dm_design()].
#' 5. Compare maps by decision quality with [dm_evaluate()] and
#'    [dm_threshold_curve()].
#' 6. Move between scales with [dm_block_matrix()], [dm_flow_matrix()],
#'    [dm_flow_from_dem()] and [dm_aggregate()].
#'
#' @section Acknowledgements:
#' The authors gratefully acknowledge the Indian Council of Agricultural
#' Research (ICAR), New Delhi, for providing the necessary research support.
#' This work was carried out as part of the institute-funded project
#' "Development of AI-driven solutions for quantifying runoff induced soil loss
#' and nutrient loss" (NRMAIISWCSIL202400300220) at ICAR-Indian Institute of
#' Soil and Water Conservation, Dehradun.
#'
#' @references
#' Raiffa, H. and Schlaifer, R. (1961). *Applied Statistical Decision Theory*.
#' Harvard University, Boston.
#'
#' Ehm, W., Gneiting, T., Jordan, A. and Krueger, F. (2016). Of quantiles and
#' expectiles: consistent scoring functions, Choquet representations and
#' forecast rankings. *Journal of the Royal Statistical Society Series B*,
#' 78(3). \doi{10.1111/rssb.12154}
#'
#' Evensen, G. (1994). Sequential data assimilation with a nonlinear
#' quasi-geostrophic model using Monte Carlo methods to forecast error
#' statistics. *Journal of Geophysical Research*, 99(C5).
#' \doi{10.1029/94JC00572}
#'
#' Meinshausen, N. (2006). Quantile regression forests. *Journal of Machine
#' Learning Research*, 7, 983-999.
#'
#' Lark, R. M., Chagumaira, C. and Milne, A. E. (2022). Decisions, uncertainty
#' and spatial information. *Spatial Statistics*, 50, 100619.
#'
#' @keywords internal
"_PACKAGE"

#' @importFrom stats pnorm dnorm qnorm rnorm runif sd var cov median quantile
#'   optim dist approx setNames cor
#' @importFrom graphics image points text par layout rect axis box legend
#'   lines abline barplot plot.new plot.window title mtext segments matplot
#' @importFrom grDevices hcl.colors adjustcolor
#' @importFrom utils head
#' @useDynLib terravai, .registration = TRUE
NULL
