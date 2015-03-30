# Run the following scripts, in order to generate all tables and figures for WCN analysis.
# Amir Shahmoradi, Sunday 11:09 PM, March 29, 2015, Wilke Lab, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_best_definition/analysis/src')

# GET INPUT DATA:
source('input_data.r')

# GENERATE QUANTILES TABLES AND THE CORRESPONDING PLOTS OF SPEARMAN RHO VS. FREE PATAMETER OF THE WCN MODEL USED:
source('get_quantiles.r')