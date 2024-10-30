# RFTMultiflagella
Model used in Lisevitch et al (2024)

This repository contains the matlab codes to replicate the computations of the model used in Lisevitch et al (2024), see preprint at: https://doi.org/10.1101/2024.03.30.587422

FlagellarNumberModel_Revision_Final.m is the main file to plot the figures

generateModelFlagellarPropulsion.m computes swimming speed and rotation rates according to parameter specifications

getBodyParamsComputeTheta.m computes the cell body friction coefficients and the modeled cell body inclination

getBodyParamsThetaCst.m computes the cell body friction coefficients for a given cell body inclination

2024-10-01_DataDFFM_Final.xlsx is the experimental data for swimming speed (v0) in µm/s, flagellar rotation rate (wf), motor rotation rate (wM) and cell body rotation rate (wb) in Hz, and their ratios, for 3 strains with different average numbers of flagella (nflag)   
