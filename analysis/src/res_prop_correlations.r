# I am writing this code to find out how strong different residue properties correlate with each other, when averaged over all PDB dataset.
# The output will be a table that contains, on each row, the mean, median and STDEV of the Spearman correlations of pairs of residue properties over the entire PDB dataset.

# Input data:
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out
#               ../../properties/res_prop_voronoiAA.out
#               ../../properties/res_prop_voronoiCA.out
#               ../../properties/res_prop_voronoiSC.out

# Last updated by Amir Shahmoradi, Tuesday 7:41 PM, Sep 30 2014, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.

res_prop_elj         = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj = res_prop_elj[!(res_prop_elj$pdb %in% excluded_pdbs),]
res_prop_elj$pdb     = factor(res_prop_elj$pdb)

res_prop_jec         = read.csv('../../jec_pdb_r4s.csv', header=T)
res_prop_jec = res_prop_jec[!(res_prop_jec$pdb %in% excluded_pdbs),]
res_prop_jec$pdb     = factor(res_prop_jec$pdb)

res_prop_hps         = read.table('../../properties/res_prop_hps.out', header=T)
res_prop_hps = res_prop_hps[!(res_prop_hps$pdb %in% excluded_pdbs),]
res_prop_hps$pdb     = factor(res_prop_hps$pdb)

res_prop_dssp        = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp = res_prop_dssp[!(res_prop_dssp$pdb %in% excluded_pdbs),]
res_prop_dssp$pdb    = factor(res_prop_dssp$pdb)

res_prop_wcn_bf      = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf = res_prop_wcn_bf[!(res_prop_wcn_bf$pdb %in% excluded_pdbs),]
res_prop_wcn_bf$pdb  = factor(res_prop_wcn_bf$pdb)

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA      = cbind(res_prop_voroAA, VAAsphericity = 4.8359758620494089221509005399179*(res_prop_voroAA$VAAvolume^(2./3.))/res_prop_voroAA$VAAarea)
res_prop_voroAA$VAAmodified_sphericity = res_prop_voroAA$VAAsphericity
res_prop_voroAA$VAAmodified_sphericity[res_prop_voroAA$VAAvolume_change_diff != 0] = -res_prop_voroAA$VAAsphericity[res_prop_voroAA$VAAvolume_change_diff != 0]
res_prop_voroAA = res_prop_voroAA[!(res_prop_voroAA$pdb %in% excluded_pdbs),]
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA      = cbind(res_prop_voroCA, VCAsphericity = 4.8359758620494089221509005399179*(res_prop_voroCA$VCAvolume^(2./3.))/res_prop_voroCA$VCAarea)
res_prop_voroCA$VCAmodified_sphericity = res_prop_voroCA$VCAsphericity
res_prop_voroCA$VCAmodified_sphericity[res_prop_voroCA$VCAvolume_change_diff != 0] = -res_prop_voroCA$VCAsphericity[res_prop_voroCA$VCAvolume_change_diff != 0]
res_prop_voroCA = res_prop_voroCA[!(res_prop_voroCA$pdb %in% excluded_pdbs),]
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC$VSCvolume^(2./3.))/res_prop_voroSC$VSCarea)
res_prop_voroSC$VSCmodified_sphericity = res_prop_voroSC$VSCsphericity
res_prop_voroSC$VSCmodified_sphericity[res_prop_voroSC$VSCvolume_change_diff != 0] = -res_prop_voroSC$VSCsphericity[res_prop_voroSC$VSCvolume_change_diff != 0]
res_prop_voroSC = res_prop_voroSC[!(res_prop_voroSC$pdb %in% excluded_pdbs),]
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)


all_res_prop = cbind( subset(res_prop_jec, select = c(pdb,r4s_JC)),
                      subset(res_prop_elj, select = c(seqent,ddgent)),
                      subset(res_prop_hps, select = c(hpskd,hpsww,hpshh)),
                      subset(res_prop_dssp, select = c(asa,rsa,hbe)),
                      subset(res_prop_wcn_bf, select = -c(pdb,resnam,resnum)),
                      subset(res_prop_voroAA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,VAAnvertices,VAAnedges,VAAvolume_change_diff,VAAvolume_change_ratio)),
                      subset(res_prop_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change_diff,VCAvolume_change_ratio)),
                      subset(res_prop_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change_diff,VSCvolume_change_ratio))
                      )


all_scors_all_pdbs = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
column_names = cbind('variable1','variable2')
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
 
  pdb_temp = all_res_prop[all_res_prop$pdb == pdb,]
  cormat = as.data.frame(cor(subset(pdb_temp, select = -c(pdb)), method='sp'))
  cormat = data.frame( variable = row.names(cormat), cormat )
  row.names(cormat) = NULL
  cormat_long = melt(cormat, id = c('variable'))
  column_name = paste0('scor_',pdb)
  column_names = cbind(column_names, column_name)
  
  if (counter == 1) { all_scors_all_pdbs = cormat_long }
  else { all_scors_all_pdbs = cbind(all_scors_all_pdbs, cormat_long$value) }
}
colnames(all_scors_all_pdbs) = column_names

write.csv( all_scors_all_pdbs, "../tables/all_scors_all_pdbs.csv", row.names=F )


# Now summarize all Spearman correlations over the entire dataset

all_scors_summary = data.frame(variable1 = all_scors_all_pdbs$variable1,
                               variable2 = all_scors_all_pdbs$variable2,
                               mean = apply(subset(all_scors_all_pdbs, select = -c(variable1,variable2)),1,mean ),
                               sd   = apply(subset(all_scors_all_pdbs, select = -c(variable1,variable2)),1,mean ),
                               as.data.frame(t(apply(subset(all_scors_all_pdbs, select = -c(variable1,variable2)),1,quantile)))
                               )

all_scors_summary$variable1 = factor(all_scors_summary$variable1)
#all_scors_summary = all_scors_summary[order(all_scors_summary[,'variable1']),]

all_scors_summary_ordered = data.frame()
for (variable1 in levels(all_scors_summary$variable1))
{
  #temp_scors = subset(all_scors_summary[all_scors_summary$variable1 == variable1, ], select=-c(variable1))
  temp_scors = all_scors_summary[all_scors_summary$variable1 == variable1, ]
  temp_scors$variable2 = factor(temp_scors$variable2)
  temp_scors = temp_scors[with(temp_scors, order(variable2)),]
  all_scors_summary_ordered = rbind(all_scors_summary_ordered, temp_scors)
}

write.csv(all_scors_summary_ordered, "../tables/all_scors_summary.csv", row.names=F )

