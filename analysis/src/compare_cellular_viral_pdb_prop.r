# I am writing this R script in an attempt to find potential structural differences between viral (ASAP) and cellular (Echave) proteins.
# Last updated by Amir Shahmoradi, Wednesday 12:23 PM, Nov 12 2014, Wilke Lab, ICMB, UT Austin

# input files:  
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_HPS_asap.out
#               ../../properties/res_prop_dssp_asap.out
#               ../../properties/res_prop_wcn_bf_asap.out


#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

pdb_prop_from_residue_prop = read.csv( "../tables/pdb_prop_from_residue_prop.csv", header=T )
pdb_prop_from_residue_prop_asap = read.csv( "../tables/pdb_prop_from_residue_prop_asap.csv", header=T )

for ()
{
    
}