# This script attempts to find out if there are any properties of proetin that could possibly affect the free parameter of WCN definitions, in particular, the power-law definition.
# Amir Shahmoradi, Thursday 9:38 PM, April 16 2015, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_best_definition/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.

pdb_CO     = read.table('../../../properties/pdb_prop_CO.out',header=T)
pdb_CO     = pdb_CO[!(pdb_CO$pdb %in% excluded_pdbs),]
pdb_CO$pdb = factor(pdb_CO$pdb)

pdb_CO_varied = read.table('../../CO_cutoff_varied/CO_cutoff_varied.out',header=F)
free_param_values = pdb_CO_varied[pdb_CO_varied$V1=='free_parameter',]
free_param_values = data.frame( t( subset(free_param_values, select=-c(V1)) ) )

pdb_CO_varied = read.table('../../CO_cutoff_varied/CO_cutoff_varied.out',header=T)
pdb_CO_varied = pdb_CO_varied[!(pdb_CO_varied$free_parameter %in% excluded_pdbs),]
pdb_CO_varied$free_parameter = factor(pdb_CO_varied$free_parameter)


sum_wcnpSC_bfSC     = read.table('../../bfactor/sum_wcnpSC_bfSC.out',header=T)
sum_wcnpSC_bfSC     = sum_wcnpSC_bfSC[!(sum_wcnpSC_bfSC$pdb %in% excluded_pdbs),]
sum_wcnpSC_bfSC$pdb = factor(sum_wcnpSC_bfSC$pdb)

sum_wcnpSC_r4sJC     = read.table('../../r4sJC/sum_wcnpSC_r4sJC.out',header=T)
sum_wcnpSC_r4sJC     = sum_wcnpSC_r4sJC[!(sum_wcnpSC_r4sJC$pdb %in% excluded_pdbs),]
sum_wcnpSC_r4sJC$pdb = factor(sum_wcnpSC_r4sJC$pdb)

sum_wcnpSC_seqent     = read.table('../../seqent/sum_wcnpSC_seqent.out',header=T)
sum_wcnpSC_seqent     = sum_wcnpSC_seqent[!(sum_wcnpSC_seqent$pdb %in% excluded_pdbs),]
sum_wcnpSC_seqent$pdb = factor(sum_wcnpSC_seqent$pdb)

sum_wcnpSC_ddgent     = read.table('../../ddgent/sum_wcnpSC_ddgent.out',header=T)
sum_wcnpSC_ddgent     = sum_wcnpSC_ddgent[!(sum_wcnpSC_ddgent$pdb %in% excluded_pdbs),]
sum_wcnpSC_ddgent$pdb = factor(sum_wcnpSC_ddgent$pdb)

cor(pdb_CO_varied$X12.000000, sum_wcnpSC_bfSC$free_param_best, method='sp' )

cor(pdb_CO$contact_orderSC, sum_wcnpSC_bfSC$free_param_best, method='sp' )
cor(pdb_CO$contact_orderSC, sum_wcnpSC_r4sJC$free_param_best, method='sp' )
cor(pdb_CO$contact_orderSC, sum_wcnpSC_seqent$free_param_best, method='sp' )
cor(pdb_CO$contact_orderSC, sum_wcnpSC_ddgent$free_param_best, method='sp' )

plot( pdb_CO$contact_orderSC, sum_wcnpSC_bfSC$free_param_best )
plot( pdb_CO$contact_orderSC, sum_wcnpSC_ddgent$free_param_best )
plot( pdb_CO_varied$X30.00000, sum_wcnpSC_ddgent$free_param_best )
plot( pdb_CO$contact_orderSC, sum_wcnpSC_r4sJC$free_param_best )
plot( pdb_CO$contact_orderSC, sum_wcnpSC_bfSC$free_param_best, ylim=c(-5,0) )

bestFP = data.frame( bfsc = sum_wcnpSC_bfSC$free_param_best
                   , r4s  = sum_wcnpSC_r4sJC$free_param_best
                   , seq  = sum_wcnpSC_seqent$free_param_best
                   , ddg  = sum_wcnpSC_ddgent$free_param_best
                   , co   = subset(pdb_CO_varied, select=-c(free_parameter))
                   )
co_cormat = as.data.frame( cor(bestFP, method='sp', ) )
#co_cormat[is.na(co_cormat)] <- 0
plot( free_param_values$X1
    , co_cormat$bfsc[-1:-4]
    , type = 'l'
    )

max(co_cormat$bfsc[-1:-4],na.rm=TRUE)

plot( free_param_values$X1
    , co_cormat$ddg[-1:-4]
    , type = 'l'
    )
plot( free_param_values$X1
    , co_cormat$r4s[-1:-4]
    , type = 'l'
    )


median(bestFP$bfsc); sd(bestFP$bfsc)
median(bestFP$r4s); sd(bestFP$r4s)
median(bestFP$seq); sd(bestFP$seq)
median(bestFP$ddg); sd(bestFP$ddg)

plot(sum_wcnpSC_bfSC$free_param_best[sum_wcnpSC_bfSC$free_param_best<5], sum_wcnpSC_ddgent$free_param_best[sum_wcnpSC_bfSC$free_param_best<5])
