# Amir Shahmoradi, Monday 3:27 PM, Oct 6 2014, Wilke Lab, ICMB, UT Austin
# This code is very similar to 'get_pdb_data.r' with the sole difference that here only select important and best representative variables are used.

#install.packages("reshape2")
library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.

pdb_CO = read.table('../../properties/pdb_prop_CO.out',header=T)
pdb_CO$pdb = factor(pdb_CO$pdb)
pdb_CO = pdb_CO[!(pdb_CO$pdb %in% excluded_pdbs),]

pdb_prop_dssp = read.table('../../properties/pdb_prop_dssp.out',header=T)
pdb_prop_dssp = cbind( pdb_prop_dssp,
                       data.frame(sum.nhbpa.dif = pdb_prop_dssp$sum.nhbps - pdb_prop_dssp$sum.nhbas),
                       data.frame(mean.nhbpa.dif = pdb_prop_dssp$mean.nhbps - pdb_prop_dssp$mean.nhbas)
                       )
pdb_prop_dssp$pdb = factor(pdb_prop_dssp$pdb)
pdb_prop_dssp = pdb_prop_dssp[!(pdb_prop_dssp$pdb %in% excluded_pdbs),]

pdb_temp = cbind( #subset(pdb_CO, select = c(pdb,natoms,contact_order,contact_orderSC,contact_orderAA)),
                  subset(pdb_CO, select = c(pdb,natoms,contact_orderSC)),
                  subset(pdb_prop_dssp, select = -c(pdb,pdb_asa))
                  )

pdb_prop_dssp_CO_long = reshape( pdb_temp,
                         idvar = 'pdb',
                         varying = names(pdb_temp[,!(names(pdb_temp) %in% 'pdb')]),
                         v.names = 'value',
                         timevar = 'variable',
                         times = names(pdb_temp[,!(names(pdb_temp) %in% 'pdb')]),
                         direction = 'long'
                         )
rownames(pdb_prop_dssp_CO_long) = NULL

# Now combine all PDB data in a single long table:

pdb_prop_from_residue_prop_select = read.csv('../tables/pdb_prop_from_residue_prop_select.csv')
all_pdb_prop_select = rbind(pdb_prop_from_residue_prop_select,pdb_prop_dssp_CO_long)

write.csv(all_pdb_prop_select, "../tables/all_pdb_prop_select.csv", row.names=F )

all_pdb_prop_select_wide = dcast(all_pdb_prop_select, pdb ~ variable, value.var = 'value', mean)
# The following lines arrange data in the chronological order of the column names.
colnames_all_pdb_prop_select_wide = data.frame( colnames = colnames( subset(all_pdb_prop_select_wide, select = -c(pdb)) ) )
colnames_all_pdb_prop_select_wide$colnames = colnames_all_pdb_prop_select_wide[with(colnames_all_pdb_prop_select_wide, order(colnames)),]
pdbs = data.frame(pdb = all_pdb_prop_select_wide$pdb)
all_pdb_prop_select_wide = cbind(pdbs, all_pdb_prop_select_wide[,as.vector(colnames_all_pdb_prop_select_wide$colnames)])
write.csv(all_pdb_prop_select_wide, "../tables/all_pdb_prop_select_wide.csv", row.names=F )

all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(pdb,sum.nssb,mean.nssb))
all_pdb_prop_cormat = as.data.frame(cor(all_pdb_prop_subset, method='spearman'))
write.csv(all_pdb_prop_cormat, "../tables/all_pdb_prop_cormat.csv", row.names=T )
write.csv(abs(all_pdb_prop_cormat), "../tables/all_pdb_prop_cormat_abs_values.csv", row.names=T )
cor(all_pdb_prop_cormat$r.rsa.seqent,all_pdb_prop_cormat$sd.vsphericitym, method='sp')
cor(all_pdb_prop_cormat$r.seqent.wcnSC,all_pdb_prop_cormat$r.seqent.vsphericitym, method='sp')
cormat_all_pdb_prop_cormat = as.data.frame(cor(all_pdb_prop_cormat, method='spearman'))
write.csv( cormat_all_pdb_prop_cormat, "../tables/cormat_all_pdb_prop_cormat.csv", row.names=T )
write.csv( abs(cormat_all_pdb_prop_cormat), "../tables/cormat_all_pdb_prop_cormat_abs_values.csv", row.names=T )




# all_pdb_prop_select$variable = factor(all_pdb_prop_select$variable)
# variable_counter = 0
# all_spcor = data.frame()
# 
# for (variable1 in levels(all_pdb_prop_select$variable))
# {
#   temp_var1 = all_pdb_prop_select[all_pdb_prop_select$variable == variable1,]
#   
#   var_spcor = data.frame()
#   start_time = proc.time()
#   for (variable2 in levels(all_pdb_prop_select$variable))
#   {
#     if (variable1 != variable2)
#     {
#       temp_var2 = all_pdb_prop_select[all_pdb_prop_select$variable == variable2,]
#       x = cor.test( temp_var1$value, temp_var2$value, method='spearman', na.action="na.omit" )
#       r = x$estimate
#       p = x$p.value
#       
#       row = data.frame( var1 = variable1,
#                         var2 = variable2,
#                         abs.spearman.r = abs(r),
#                         spearman.r = r,
#                         spearman.p = p
#                         )
#       
#       var_spcor = rbind(var_spcor,row)
#     }
#   }
#   
#   var_spcor_ordered = var_spcor[with(var_spcor, order(-abs.spearman.r)),]
#   
#   filename = paste0( '../tables/correlations/select_variables/',variable1,'.csv' )
#   write.csv( var_spcor_ordered, file = filename, row.names=F )
#   
#   all_spcor = rbind (all_spcor, var_spcor_ordered)
#   
#   variable_counter = variable_counter + 1
#   cat(variable_counter, variable1, nrow(temp_var1), 'time taken: ', proc.time()-start_time, '\n')
#   
# }
# 
# write.csv( all_spcor, file = '../tables/correlations/select_variables/all_correlations.csv', row.names=F )
