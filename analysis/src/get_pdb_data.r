# Amir Shahmoradi, Thursday 4:20 PM, Aug 21 2014, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

pdb_CO = read.table('../../properties/pdb_prop_CO.out',header=T)

cor.test(all_pdb_prop$r.seqent_wcnca,pdb_CO$contact_orderSC, method='spearman')
cor.test(all_pdb_prop$r.ddgent_wcnca,pdb_CO$contact_order, method='spearman')
cor.test(all_pdb_prop$mwcn,pdb_CO$contact_order, method='spearman')
cor.test(all_pdb_prop$nhbon,pdb_CO$contact_orderSC, method='spearman')
cor.test(all_pdb_prop$nhbon,all_pdb_prop$mrsa, method='spearman')
cor.test(all_pdb_prop$nhbon,all_pdb_prop$mwcn, method='spearman')
cor.test(all_pdb_prop$nhbas,all_pdb_prop$r.seqent_wcnca, method='spearman')
cor.test(all_pdb_prop$nhbps,all_pdb_prop$r.seqent_rsa, method='spearman')
cor.test(all_pdb_prop$r.ddgent_wcnca,all_pdb_prop$vwcn, method='spearman')
cor.test(pdb_CO$contact_orderAA,pdb_CO$contact_orderSC, method='spearman')
