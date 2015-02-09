# This is the main R script that runs all R scripts in the proper order, to collect, generate, and analyze all data and output results, tables and figures:
# Last updated by Amir Shahmoradi, Sunday 7:54 PM, Feb 8 2015, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# First collect all residue-level data, that is, all site specific properties.
source('get_res_data.r')

# Now make plots of correlation of Voronoi cell area with other site-specific variables for all proteins in the form of box plots.
 
