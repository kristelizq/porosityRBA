# porosityRBA
Residual Bouguer Anomaly (RBA) data and functions to plot best fit models of the data and corresponding Bayesian Information Criteria (BIC) of each model.

Data:
RBArobbins18Goossens20_rho2550_taper.csv contains the RBA data used in Izquierdo et al (2021) for lunar craters with diameters between 10 and 30 km.
regionCraters.csv contains the location, size and region of the lunar craters used to obtain the RBA. The regions are SPA, highlands or mare.

Functions:
BICmodels.m and BICmodels_regions.m calculate the BIC of one and two slopes models for the whole data RBA data or for the RBA data per region.
rba2s2fModel.m computes the best-fit two slopes models for the whole RBA data set with uncertainty.
