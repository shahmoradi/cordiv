# This R code performs regularized regression on pdb properties data.
# Amir Shahmoradi, Sunday 4:47 AM, Oct 19 2014, Wilke Lab, ICMB, UT Austin

install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(pdb,sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
seqent_scors = grep('^r.+seqent',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
r4sJC_scors = grep('^r.+r4sJC',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
not_in_pcr = c(seqent_scors,r4sJC_scors)

fit = glmnet(x = as.matrix(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)]), y = as.vector(all_pdb_prop_subset$r.seqent.wcnSC))
predictor_names = colnames(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)])
predictor_names[145]
plot(fit, label = TRUE)
print(fit)
coef(fit,s=0.1)

