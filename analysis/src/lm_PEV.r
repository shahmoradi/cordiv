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



# Now generate a single file containing all coefficients for the 4 most important structure-sequence correlations:

# Fisrt seqent-structure correlations
all_ssscors = abs(subset(all_pdb_prop_cormat, select = names(all_pdb_prop_cormat) %in% cornames_seqent))
all_ssscors = data.frame(variable = all_pdb_prop_cormat$X,all_ssscors)
all_ssscors = all_ssscors[with(all_ssscors, order(-r.seqent.wcnSC)),]
#all_ssscors = all_ssscors[with(all_ssscors, order(-r.rsa.seqent)),]

#plot( log10(order(-all_ssscors$r.seqent.wcnSC))
pdf( "../figures/PEV_structure_seqent_cors_rank_data_wcnSC_ordered.pdf", width=4.5, height=4, useDingbats=FALSE )
#pdf( "../figures/PEV_structure_seqent_cors_rank_data_RSA_ordered.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(
    log10(order(-all_ssscors$r.seqent.wcnSC))
    #log10(order(-all_ssscors$r.rsa.seqent))
    , all_ssscors$r.seqent.wcnSC^2
    #, all_ssscors$r.rsa.seqent^2
    , type = 'l'
    #, xlim = c(0,3)
    , ylim = c(0,0.3)
    , xlab = "Log10 (Order of Explanatory Variables)"
    , ylab = "Proportion of Explained Variance (Rank Data)"
    )
lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.seqent.varea^2, col = 'blue')
#lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.seqent.varea^2, col = 'blue')
lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.rsa.seqent^2, col = 'red')
#lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.seqent.wcnSC^2, col = 'red')
lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.ddgent.seqent^2, col = 'green')
#lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.ddgent.seqent^2, col = 'green')
##lines(log10(order(-all_ssscors$r.seqent.wcnSC)), all_ssscors$r.bfSC.seqent, col = 'orange')
##lines(log10(order(-all_ssscors$r.rsa.seqent)), all_ssscors$r.bfSC.seqent^2, col = 'orange')

legend( 0.0, 0.32
        , c('r.seqent.wcnSC','r.seqent.varea')
        #, c('r.seqent.rsa','r.seqent.varea')
        #, pch=19
        , lty=c(1,1)
        , col=c('black','blue')
        , bty='n', cex=0.9)
legend( 1.0, 0.32
        , c('r.seqent.rsa', 'r.seqent.ddgent') #, 'r.seqent.bfSC')
        #, c('r.seqent.wcnSC', 'r.seqent.ddgent') #, 'r.seqent.bfSC')
        #, pch=19
        , lty=c(1,1)
        , col=c('red','green')  #,'orange')
        , bty='n', cex=0.9)

dev.off()



# Then r4sJC-structure correlations
all_ssscors = abs(subset(all_pdb_prop_cormat, select = names(all_pdb_prop_cormat) %in% cornames_r4sJC))
all_ssscors = data.frame(variable = all_pdb_prop_cormat$X,all_ssscors)
#all_ssscors = all_ssscors[with(all_ssscors, order(-r.seqent.wcnSC)),]
all_ssscors = all_ssscors[with(all_ssscors, order(-r.r4sJC.rsa)),]

#plot( log10(order(-all_ssscors$r.seqent.wcnSC))
pdf( "../figures/PEV_structure_r4sJC_cors_rank_data.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot( log10(order(-all_ssscors$r.r4sJC.rsa))
      , all_ssscors$r.r4sJC.rsa^2
      , type = 'l'
      #, xlim = c(0,3)
      , ylim = c(0,0.15)
      , xlab = "Log10 (Order of Explanatory Variables)"
      , ylab = "Proportion of Explained Variance (Rank Data)"
)
#lines(log10(order(-all_ssscors$r.r4sJC.wcnSC)), all_ssscors$r.r4sJC.varea, col = 'blue')
lines(log10(order(-all_ssscors$r.r4sJC.rsa)), all_ssscors$r.r4sJC.varea^2, col = 'blue')
#lines(log10(order(-all_ssscors$r.r4sJC.wcnSC)), all_ssscors$r.r4sJC.wcnSC, col = 'red')
lines(log10(order(-all_ssscors$r.r4sJC.rsa)), all_ssscors$r.r4sJC.wcnSC^2, col = 'red')
#lines(log10(order(-all_ssscors$r.r4sJC.wcnSC)), all_ssscors$r.ddgent.r4sJC, col = 'green')
lines(log10(order(-all_ssscors$r.r4sJC.rsa)), all_ssscors$r.ddgent.r4sJC^2, col = 'green')
#lines(log10(order(-all_ssscors$r.r4sJC.wcnSC)), all_ssscors$r.bfSC.r4sJC, col = 'orange')
#lines(log10(order(-all_ssscors$r.r4sJC.rsa)), all_ssscors$r.bfSC.r4sJC^2, col = 'orange')

legend( 0.0, 0.32, c('r.r4sJC.rsa','r.r4sJC.varea')
        #, pch=19
        , lty=c(1,1)
        , col=c('black','blue')
        , bty='n', cex=0.9)
legend( 1.0, 0.32, c('r.r4sJC.wcnSC', 'r.r4sJC.ddgent') #, 'r.r4sJC.bfSC')
        #, pch=19
        , lty=c(1,1)
        , col=c('red','green')  #,'orange')
        , bty='n', cex=0.9)

dev.off()



