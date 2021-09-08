# porosityRBA
Residual Bouguer Anomaly (RBA) data and functions to plot best fit models of the data and corresponding Bayesian Information Criteria (BIC) of each model.

Data:

RBArobbins18Goossens20_rho2550_taper.csv contains the RBA data used in Izquierdo et al (submitted) for lunar craters with diameters between 10 and 30 km.

regionCraters.csv contains the location, size and region of the lunar craters between 10 and 30 km (Robbins,2019). The regions are SPA, highlands or mare.

Functions:

BICmodels.m and BICmodels_regions.m calculate the BIC of one and two slopes models for the whole RBA data or for the RBA data per region. No inputs.

rba2s2fModel.m computes the best-fit two slopes models for the whole RBA data set with uncertainty. Inputs: number of iterations for bootstrapping and option of printing figures 'yes' or 'nop'.

bestModelUnc_regions.m computes the best-fit model for the RBA of each region with uncertainty. Inputs: number of iterations for bootstrapping and option of printing figures 'yes' or 'nop'.
