# This code generates a ranked version of the dataframe containing all residue PDB data.
# Amir Shahmoradi, Thursday 6:26 PM, Feb 26 2014, Wilke Lab, ICMB, UT Austin 

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

pdb_temp = data.frame()
best_structural_predictors_of_ER = data.frame()  # This data frame will contain correlations of select structural variables with r4sJC evolutionary rates, for each pdb file on eaxch row
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(pdb,seqent,ddgent))
                    , subset(pdb_jec, select = c(zr4s_JC))
                    , subset(pdb_hps, select = c(hpshh))
                    , subset(pdb_dssp, select = c(rsa,hbe))
                    , subset(pdb_wcn_bf, select = c(wcnSC,wcnCA,bfSC))
                    , subset(pdb_voroSC, select = c(VSCarea))
                    , subset(pdb_voroCA, select = c(VCAarea))
  )
  r.rsa.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$rsa,method='sp')
  r.wcnSC.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$wcnSC,method='sp')
  r.wcnCA.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$wcnCA,method='sp')
  r.vareaSC.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$VSCarea,method='sp')
  r.ddgent.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$ddgent,method='sp')
  r.vareaCA.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$VCAarea,method='sp')
  r.bfSC.r4sJC = cor(pdb_temp$zr4s_JC,pdb_temp$bfSC,method='sp')
  
  #row = data.frame( pdb = pdb, rsa = r.rsa.r4sJC, wcn = r.wcn.r4sJC, vareaSC = r.vareaSC.r4sJC, ddgent = r.ddgent.r4sJC )
  row = data.frame( pdb, r.wcnCA.r4sJC, r.wcnSC.r4sJC, r.rsa.r4sJC, r.vareaSC.r4sJC, r.ddgent.r4sJC , r.vareaCA.r4sJC , r.bfSC.r4sJC )
  best_structural_predictors_of_ER = rbind( best_structural_predictors_of_ER, row )
}

write.csv(best_voronoi_predictors_of_ER, file = "../tables/best_structural_predictors_of_ER.csv", row.names=F )












all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.
all_pdb_prop_select_wide_rank = all_pdb_prop_subset$pdb

for (column in colnames(all_pdb_prop_subset))
{
  if (column != 'pdb')
  {
    temp_column = data.frame(rank(all_pdb_prop_subset[[column]]))
    all_pdb_prop_select_wide_rank = cbind(all_pdb_prop_select_wide_rank, temp_column)
  }
}

colnames(all_pdb_prop_select_wide_rank) = colnames(all_pdb_prop_subset)

write.csv(all_pdb_prop_select_wide_rank, "../tables/all_pdb_prop_select_wide_rank.csv", row.names=F )

all_pdb_prop_rank_cormat = as.data.frame(cor(subset(all_pdb_prop_select_wide_rank, select=-c(pdb))), method='spearman')
write.csv(all_pdb_prop_rank_cormat, "../tables/all_pdb_prop_rank_cormat.csv", row.names=T )

all_pdb_prop_select_wide_rank = read.csv("../tables/all_pdb_prop_select_wide_rank.csv", header=T)

#regressors = subset(all_pdb_prop_select_wide_rank, select=-c(pdb,r.seqent.wcnSC))
#pcrdata = cbind(data.frame( r.seqent.wcnSC = all_pdb_prop_select_wide_rank$r.seqent.wcnSC), regressors)
#pcrdata_names = colnames(pcrdata)

#pcr_out = pcr(r.seqent.wcnSC ~ . , data = pcrdata, y = T)
#summary(pcr_out)

#summary(lm (r.seqent.wcnSC ~ . , all_pdb_prop_subset))
#summary(lm (r.seqent.wcnSC ~ . , all_pdb_prop_rank_cormat))
#pca_out = prcomp(pcr_data)
#write.csv(pca_out$, '../tables/pca_loads.csv', row.names = T)
