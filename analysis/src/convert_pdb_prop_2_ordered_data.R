# This code generates an ordered version of the dataframe containing all PDB data.
# Amir Shahmoradi, Tuesday 1:27 AM, Oct 7 2014, Wilke Lab, ICMB, UT Austin 

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
all_pdb_prop_select_wide_order = all_pdb_prop_subset$pdb

for (column in colnames(all_pdb_prop_subset))
{
  if (column != 'pdb')
  {
    temp_column = data.frame(rank(all_pdb_prop_subset[[column]]))
    all_pdb_prop_select_wide_order = cbind(all_pdb_prop_select_wide_order, temp_column)
  }
}

colnames(all_pdb_prop_select_wide_order) = colnames(all_pdb_prop_subset)

write.csv(all_pdb_prop_select_wide_order, "../tables/all_pdb_prop_select_wide_order.csv", row.names=F )

all_pdb_prop_order_cormat = as.data.frame(cor(subset(all_pdb_prop_select_wide_order, select=-c(pdb))), method='spearman')
write.csv(all_pdb_prop_order_cormat, "../tables/all_pdb_prop_order_cormat.csv", row.names=T )

regressors = subset(all_pdb_prop_select_wide_order, select=-c(pdb,r.seqent.wcnSC))
pcrdata = cbind(data.frame( r.seqent.wcnSC = all_pdb_prop_select_wide_order$r.seqent.wcnSC), regressors)
pcrdata_names = colnames(pcrdata)

pcr_out = pcr(r.seqent.wcnSC ~ . , data = pcrdata, y = T)
summary(pcr_out)

summary(lm (r.seqent.wcnSC ~ . , all_pdb_prop_subset))
summary(lm (r.seqent.wcnSC ~ . , all_pdb_prop_order_cormat))
pca_out = prcomp(pcr_data)
write.csv(pca_out$, '../tables/pca_loads.csv', row.names = T)
