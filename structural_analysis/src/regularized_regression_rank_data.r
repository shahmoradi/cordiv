# This R code performs regularized regression on pdb properties RANK data.
# Amir Shahmoradi, Friday 8:09 PM, Oct 24 2014, Wilke Lab, ICMB, UT Austin

# install.packages("glmnet", repos = "http://cran.us.r-project.org")
library(glmnet)

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select_wide_rank = read.csv("../tables/all_pdb_prop_select_wide_rank.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide_rank, select = -c(pdb)) # ,sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
seqent_scors = grep('^r.+seqent',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
r4sJC_scors = grep('^r.+r4sJC',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
not_in_pcr = c(seqent_scors,r4sJC_scors)

x = as.matrix(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)])   # predictor variables

for (correlation in not_in_pcr)
{
  cat (paste(' Fitting for correlation:', correlation, '\n'))
  fit = glmnet(x = x,
               y = as.vector(all_pdb_prop_subset[[correlation]]),
               alpha = 0.5   # elastic-net mixing parameter, with a range from 0 (ridge regression) to 1 (lasso regression).
               )

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
  
  filename = paste0('../tables/regularized_regression_rank_data/',correlation,'_lmin_rank.csv')
  write.csv(lambda.min.coef, filename, row.names = F)
  
  filename = paste0('../tables/regularized_regression_rank_data/',correlation,'_l1se_rank.csv')
  write.csv(lambda.1se.coef, filename, row.names = F)
  
  # Now generate the error-lambda plots:
  
  # filename = paste0('../figures/regularized_regression_rank_data/',correlation,'.pdf')
  # pdf( filename, width=6, height=5, useDingbats=FALSE )
  # plot(cvfit)
  # dev.off()  
  # cvfit$lambda.min
  #colnames(all_pdb_prop_subset[177])
  #plot(fit, xvar = "dev", label = TRUE)
}

# Now generate a single file containing all coefficients for the 5 most important structure-sequence correlations:

cornames_r4sJC  = c('r.bfSC.r4sJC','r.ddgent.r4sJC','r.r4sJC.rsa','r.r4sJC.varea','r.r4sJC.wcnSC')

filename = paste0('../tables/regularized_regression_rank_data/r.r4sJC.wcnSC_lmin_rank.csv')
r.wcnSC_lmin_rank = read.csv(filename, header = T)
r.wcnSC_lmin_rank = r.wcnSC_lmin_rank[r.wcnSC_lmin_rank$variable != '(Intercept)',]
r.wcnSC_lmin_rank = r.wcnSC_lmin_rank[with(r.wcnSC_lmin_rank, order(variable)),]

filename = paste0('../tables/regularized_regression_rank_data/r.r4sJC.varea_lmin_rank.csv')
r.varea_lmin_rank = read.csv(filename, header = T)
r.varea_lmin_rank = r.varea_lmin_rank[r.varea_lmin_rank$variable != '(Intercept)',]
r.varea_lmin_rank = r.varea_lmin_rank[with(r.varea_lmin_rank, order(variable)),]

filename = paste0('../tables/regularized_regression_rank_data/r.r4sJC.rsa_lmin_rank.csv')
r.rsa_lmin_rank = read.csv(filename, header = T)
r.rsa_lmin_rank = r.rsa_lmin_rank[r.rsa_lmin_rank$variable != '(Intercept)',]
r.rsa_lmin_rank = r.rsa_lmin_rank[with(r.rsa_lmin_rank, order(variable)),]

filename = paste0('../tables/regularized_regression_rank_data/r.ddgent.r4sJC_lmin_rank.csv')
r.ddgent_lmin_rank = read.csv(filename, header = T)
r.ddgent_lmin_rank = r.ddgent_lmin_rank[r.ddgent_lmin_rank$variable != '(Intercept)',]
r.ddgent_lmin_rank = r.ddgent_lmin_rank[with(r.ddgent_lmin_rank, order(variable)),]

filename = paste0('../tables/regularized_regression_rank_data/r.bfSC.r4sJC_lmin_rank.csv')
r.bfSC_lmin_rank = read.csv(filename, header = T)
r.bfSC_lmin_rank = r.bfSC_lmin_rank[r.bfSC_lmin_rank$variable != '(Intercept)',]
r.bfSC_lmin_rank = r.bfSC_lmin_rank[with(r.bfSC_lmin_rank, order(variable)),]

all_lmin_rank_coef = data.frame( variable    = r.wcnSC_lmin_rank$variable
                               , wcnSC.coef  = r.wcnSC_lmin_rank$abs.coef
                               , varea.coef  = r.varea_lmin_rank$abs.coef
                               , rsa.coef    = r.rsa_lmin_rank$abs.coef
                               , ddgent.coef = r.ddgent_lmin_rank$abs.coef
                               , bfSC.coef   = r.bfSC_lmin_rank$abs.coef
                               )
all_lmin_rank_coef = all_lmin_rank_coef[with(all_lmin_rank_coef, order(-wcnSC.coef)),]

plot(log10(order(-all_lmin_rank_coef$wcnSC.coef)), all_lmin_rank_coef$wcnSC.coef)
points(log10(order(-all_lmin_rank_coef$wcnSC.coef)), all_lmin_rank_coef$varea.coef, col = 'red')
points(log10(order(-all_lmin_rank_coef$wcnSC.coef)), all_lmin_rank_coef$rsa.coef, col = 'blue')





# plot(all_lmin_rank_coef$wcnSC.coef, all_lmin_rank_coef$rsa.coef)
# plot(all_lmin_rank_coef$wcnSC.coef, all_lmin_rank_coef$varea.coef)
# 
# 
# filename = paste0('../tables/regularized_regression_rank_data/',cornames_r4sJC[5],'_lmin_rank.csv')
# r.lmin_rank = read.csv(filename, header = T)
# r.lmin_rank = r.lmin_rank[r.lmin_rank$variable != '(Intercept)',]
# #r.r4sJC_wcnSC_lmin_rank = r.r4sJC_wcnSC_lmin_rank[r.r4sJC_wcnSC_lmin_rank$coef != 0.,]
# plot(log10(order(-r.lmin_rank$abs.coef)), r.lmin_rank$abs.coef)







# predictor_names = colnames(all_pdb_prop_subset[,-which(names(all_pdb_prop_subset) %in% not_in_pcr)])
# predictor_names[103]
# plot(fit, label = TRUE)
# plot(fit, xvar = 'dev', label = TRUE,
#      xlim = c(0.,0.1),
#      ylim = c(-0.05,0.05)
#      )
# print(fit)
# coef(fit,s=0.1)

