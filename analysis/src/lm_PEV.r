# This R script performs linear regression of sequence-structure on the set of variables that have the strongest explanatory power according to regularized regression results.
# Last updated by Amir Shahmoradi, Sunday 6:33 PM, November 2 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(pdb,sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
all_pdb_prop_cormat = read.csv("../tables/all_pdb_prop_cormat.csv", header=T )
all_pdb_prop_cormat_pvalues = read.csv("../tables/all_pdb_prop_cormat_pvalues.csv", header=T )
#plot(log10(all_pdb_prop_cormat_pvalues$contact_orderSC),abs(all_pdb_prop_cormat$contact_orderSC))

seqent_scors = grep('^r.+seqent',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
r4sJC_scors = grep('^r.+r4sJC',colnames(all_pdb_prop_subset),value = TRUE)  # Select column names that correspond to seqent-variable correlations
ssscors = c(seqent_scors,r4sJC_scors)   # structure-sequence spearman correlation names
'%ni%' = Negate('%in%')  # not in

#all_pdb_prop_cormat = subset(all_pdb_prop_cormat, select = names(all_pdb_prop_cormat) %ni% ssscors)
all_pdb_prop_cormat = all_pdb_prop_cormat[!(all_pdb_prop_cormat$X %in% ssscors),]

#all_pdb_prop_cormat_pvalues = subset(all_pdb_prop_cormat_pvalues, select = names(all_pdb_prop_cormat_pvalues) %ni% ssscors)
all_pdb_prop_cormat_pvalues = all_pdb_prop_cormat_pvalues[!(all_pdb_prop_cormat_pvalues$X %in% ssscors),]

for (correlation in ssscors)
{
  nvar = length(all_pdb_prop_cormat[[correlation]][all_pdb_prop_cormat_pvalues[[correlation]] <= 0.01])
  cat (correlation, nvar, '\n')
}

# Now find the PEV by the variance of sequence entropy and r4s:

cornames_seqent = c('r.bfSC.seqent','r.ddgent.seqent','r.rsa.seqent','r.seqent.varea','r.seqent.wcnSC')
cornames_r4sJC  = c('r.bfSC.r4sJC','r.ddgent.r4sJC','r.r4sJC.rsa','r.r4sJC.varea','r.r4sJC.wcnSC')

regressor = subset(all_pdb_prop_subset, select = colnames(all_pdb_prop_subset)[all_pdb_prop_cormat_pvalues[[cornames_seqent[5]]] <= 0.01])

lfit = lm(all_pdb_prop_subset$r.seqent.wcnSC ~ ., regressor)



# Now generate a single file containing all coefficients for the 5 most important structure-sequence correlations:

all_ssscors = abs(subset(all_pdb_prop_cormat, select = names(all_pdb_prop_cormat) %in% cornames_seqent))
all_ssscors = data.frame(variable = all_pdb_prop_cormat$X,all_ssscors)
#all_ssscors = all_ssscors[with(all_ssscors, order(-r.seqent.wcnSC)),]
all_ssscors = all_ssscors[with(all_ssscors, order(-r.rsa.seqent)),]

#plot( log10(order(-all_ssscors$r.seqent.wcnSC))
plot( log10(order(-all_ssscors$r.rsa.seqent))
    , all_ssscors$r.rsa.seqent
    , type = 'l'
    #, xlim = c(0,3)
    , ylim = c(0,0.50)
    )
#lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.seqent.varea, col = 'blue')
lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.seqent.varea, col = 'blue')
#lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.seqent.wcnSC, col = 'red')
lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.seqent.wcnSC, col = 'red')
#lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.ddgent.seqent, col = 'green')
lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.ddgent.seqent, col = 'green')
#lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.bfSC.seqent, col = 'orange')
lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.bfSC.seqent, col = 'orange')

plot(order(-all_ssscors$r.seqent.wcnSC), all_ssscors$r.seqent.wcnSC, type = 'l')
lines(order(-all_ssscors$r.seqent.wcnSC), all_ssscors$r.seqent.varea, col = 'blue')
lines(order(-all_ssscors$r.seqent.wcnSC), all_ssscors$r.rsa.seqent, col = 'red')
lines(order(-all_ssscors$r.seqent.wcnSC), all_ssscors$r.ddgent.seqent, col = 'green')


