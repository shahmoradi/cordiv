# This R code performs regularized regression on pdb properties data.
# Amir Shahmoradi, Sunday 4:47 AM, Oct 19 2014, Wilke Lab, ICMB, UT Austin

# install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(pdb,sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
seqent_scors = grep('^r.+seqent',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
r4sJC_scors = grep('^r.+r4sJC',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
not_in_pcr = c(seqent_scors,r4sJC_scors)

x = as.matrix(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)])   # predictor variables

for (correlation in not_in_pcr)
{
  cat (paste(' Now fitting for correlation:', correlation, '\n'))
  #fit = glmnet(x = x,
  #             y = as.vector(all_pdb_prop_subset[[correlation]]),
  #             alpha = 0.5   # elastic-net mixing parameter, with a range from 0 (ridge regression) to 1 (lasso regression).
  #             )

  cvfit = cv.glmnet(x = x,
                    y = as.vector(all_pdb_prop_subset[[correlation]]),
                    alpha = 0.5   # elastic-net mixing parameter, with a range from 0 (ridge regression) to 1 (lasso regression).
                    )
  
  lambda.min.coef = data.frame(as.matrix(coef(cvfit, s = 'lambda.min')))
  lambda.min.coef = data.frame(variable = row.names(lambda.min.coef),
                               coef     = lambda.min.coef$X1,
                               abs.coef = abs(lambda.min.coef$X1)
                               )
  lambda.min.coef = lambda.min.coef[with(lambda.min.coef, order(-abs.coef)),]

  lambda.1se.coef = data.frame(as.matrix(coef(cvfit, s = 'lambda.1se')))
  lambda.1se.coef = data.frame(variable = row.names(lambda.1se.coef),
                               coef     = lambda.1se.coef$X1,
                               abs.coef = abs(lambda.1se.coef$X1)
                               )
  lambda.1se.coef = lambda.1se.coef[with(lambda.1se.coef, order(-abs.coef)),]
  
  filename = paste0('../tables/regularized_regression/',correlation,'_lmin.csv')
  write.csv(lambda.min.coef, filename, row.names = F)
  
  filename = paste0('../tables/regularized_regression/',correlation,'_l1se.csv')
  write.csv(lambda.1se.coef, filename, row.names = F)
  
  # Now generate the error-lambda plots:
  
  filename = paste0('../figures/regularized_regression/',correlation,'.pdf')
  pdf( filename, width=6, height=5, useDingbats=FALSE )
  plot(cvfit)
  dev.off()  
  #cvfit$lambda.min
}



predictor_names = colnames(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)])
predictor_names[103]
plot(fit, label = TRUE)
plot(fit, xvar = 'lambda', label = TRUE)
print(fit)
coef(fit,s=0.1)

