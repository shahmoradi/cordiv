# This is the main R script that runs all R scripts in the proper order, to collect, generate, and analyze all data and output results, tables and figures:
# Last updated by Amir Shahmoradi, Sunday 7:54 PM, Feb 8 2015, Wilke Lab, ICMB, UT Austin

# install.packages("reshape2")
# install.packages("corrplot")
# install.packages('deldir')  # used by voronoi_diagram.r
# install.packages('Rpdb')  # used by getcrd.r in voronoi_diagram.r
# library("reshape2")
# library('corrplot')
# library('deldir')
# library('Rpdb')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

npdbs = 209 # number of pdbs in dataset. Will be used by other subroutines below

# First collect all residue-level data, that is, all site specific properties.
source('get_res_data.r')

# Generate data for the best definition of Voronoi cell areas.
# Make plots of correlation of Voronoi cell area with other site-specific variables for all proteins in the form of box plots.
source('best_varea_select_variables.r')

# Generate data for the best definition of B factors:
# Make plots of correlation of B factor with other site-specific variables for all proteins in the form of box plots.
source('best_bfactor_select_variables.r')

# Generate data for the best definition of WCN:
# Make plots of correlation of WCN with other site-specific variables for all proteins in the form of box plots.
source('best_wcn_select_variables.r')

# Generate an example Voronoi diagram
source('voronoi_diagram.r')

# Generate density plots comparing the predictive power of Voronoi cell area, Weighted Contact Number, and RSA about sequence evolution:
source('best_structural_predictors_of_ER.r')    # for Evolutionary Rates: r4sJC
source('best_structural_predictors_of_SE.r')    # for Sequence Entropy: seqent

# Generate density plots comparing the predictive power of different properties of Voronoi cell:
source('best_voronoi_predictors_of_ER.r')

 
