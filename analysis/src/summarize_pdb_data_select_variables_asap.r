# Amir Shahmoradi, Wednesday 3:31 PM, Nov 12 2014, Wilke Lab, ICMB, UT Austin
# This code is very similar to 'get_pdb_data.r' with the sole difference that here only select important and best representative variables are used.

#install.packages("reshape2")
library("reshape2")
#install.packages("Hmisc")
#library("Hmisc")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

nonviral_pdbs = c('1AJ8_A','1AOR_A','1CTS_A','1MP9_A','3GSZ_A','3I5K_A')   # These are thermophilic proteins, in addition to the 2 viral PDBs that are from the same family as 3GOL_A and therefore redundant.

# pdb_general = read.csv('../../properties/pdb_prop_general.csv',header=T)
# pdb_general = data.frame(pdb  = paste0(pdb_general$Protein,'_',pdb_general$chain),
#                          sres = pdb_general$Resolution,    # structure resolution
#                          nseq = pdb_general$Number.of.Sequences   # number of sequences in the alignment
#                          )
# pdb_general = pdb_general[!(pdb_general$pdb %in% nonviral_pdbs),]
# pdb_general$pdb = factor(pdb_general$pdb)

pdb_CO_asap = read.table('../../properties/pdb_prop_CO_asap.out',header=T)
pdb_CO_asap = pdb_CO_asap[!(pdb_CO_asap$pdb %in% nonviral_pdbs),]

pdb_prop_dssp_asap = read.table('../../properties/pdb_prop_dssp_asap.out',header=T)
pdb_prop_dssp_asap = cbind( pdb_prop_dssp_asap,
                       data.frame(sum.nhbpa.dif = pdb_prop_dssp_asap$sum.nhbps - pdb_prop_dssp_asap$sum.nhbas),
                       data.frame(mean.nhbpa.dif = pdb_prop_dssp_asap$mean.nhbps - pdb_prop_dssp_asap$mean.nhbas)
                       )
pdb_prop_dssp_asap = pdb_prop_dssp_asap[!(pdb_prop_dssp_asap$pdb %in% nonviral_pdbs),]
pdb_prop_dssp_asap$pdb = factor(pdb_prop_dssp_asap$pdb)

pdb_temp = cbind( #subset(pdb_CO_asap, select = c(pdb,natoms,contact_order,contact_orderSC,contact_orderAA)),
                  #pdb_general,
                  subset(pdb_CO_asap, select = c(pdb,natoms,contact_orderSC)),
                  subset(pdb_prop_dssp_asap, select = -c(pdb,pdb_asa))
                  )

pdb_prop_dssp_asap_CO_long = reshape( pdb_temp,
                         idvar = 'pdb',
                         varying = names(pdb_temp[,!(names(pdb_temp) %in% 'pdb')]),
                         v.names = 'value',
                         timevar = 'variable',
                         times = names(pdb_temp[,!(names(pdb_temp) %in% 'pdb')]),
                         direction = 'long'
                         )
rownames(pdb_prop_dssp_asap_CO_long) = NULL

# Now combine all PDB data in a single long table:

pdb_prop_from_residue_prop_select_asap = read.csv('../tables/pdb_prop_from_residue_prop_select_asap.csv', header = T)
all_pdb_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,pdb_prop_dssp_asap_CO_long)

write.csv(all_pdb_prop_select_asap, "../tables/all_pdb_prop_select_asap.csv", row.names=F )

all_pdb_prop_select_wide_asap = dcast(all_pdb_prop_select_asap, pdb ~ variable, value.var = 'value')
# The following lines arrange data in the chronological order of the column names.
colnames_all_pdb_prop_select_wide = data.frame( colnames = colnames( subset(all_pdb_prop_select_wide_asap, select = -c(pdb)) ) )
colnames_all_pdb_prop_select_wide$colnames = colnames_all_pdb_prop_select_wide[with(colnames_all_pdb_prop_select_wide, order(colnames)),]
pdbs = data.frame(pdb = all_pdb_prop_select_wide_asap$pdb)
all_pdb_prop_select_wide_asap = cbind(pdbs, all_pdb_prop_select_wide_asap[,as.vector(colnames_all_pdb_prop_select_wide$colnames)])
write.csv(all_pdb_prop_select_wide_asap, "../tables/all_pdb_prop_select_wide_asap.csv", row.names=F )
