# Calculates the spearman correlation between Sequence Entropy (seqent) and the MODIFIED Voronoi volumes of the Amino Acid.

# Amir Shahmoradi, Friday 5:58 PM, Sep 19 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

#res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
#res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)
#res_prop_voroAA      = cbind(res_prop_voroAA, VAAmodified_volume_diff = res_prop_voroAA$VAAvolume, VAAmodified_volume_ratio = res_prop_voroAA$VAAvolume)
#maxval = max(res_prop_voroAA$VAAvolume)
#res_prop_voroAA$VAAmodified_volume_diff[res_prop_voroAA$VAAvolume_change_diff != 0] = maxval
#res_prop_voroAA$VAAmodified_volume_diff  = res_prop_voroAA$VAAmodified_volume_diff + res_prop_voroAA$VAAvolume_change_diff
#res_prop_voroAA$VAAmodified_volume_ratio[res_prop_voroAA$VAAvolume_change_ratio != 1] = maxval
#res_prop_voroAA$VAAmodified_volume_ratio = res_prop_voroAA$VAAmodified_volume_ratio * res_prop_voroAA$VAAvolume_change_ratio
#
#res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
#res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)
#res_prop_voroCA      = cbind(res_prop_voroCA, VCAmodified_volume_diff = res_prop_voroCA$VCAvolume, VCAmodified_volume_ratio = res_prop_voroCA$VCAvolume)
#maxval = max(res_prop_voroCA$VCAvolume)
#res_prop_voroCA$VCAmodified_volume_diff[res_prop_voroCA$VCAvolume_change_diff != 0] = maxval
#res_prop_voroCA$VCAmodified_volume_diff  = res_prop_voroCA$VCAmodified_volume_diff + res_prop_voroCA$VCAvolume_change_diff
#res_prop_voroCA$VCAmodified_volume_ratio[res_prop_voroCA$VCAvolume_change_ratio != 1] = maxval
#res_prop_voroCA$VCAmodified_volume_ratio = res_prop_voroCA$VCAmodified_volume_ratio * res_prop_voroCA$VCAvolume_change_ratio

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCmodified_volume_diff = res_prop_voroSC$VSCvolume, VSCmodified_volume_ratio = res_prop_voroSC$VSCvolume)
maxval = max(res_prop_voroSC$VSCvolume)
res_prop_voroSC$VSCmodified_volume_diff[res_prop_voroSC$VSCvolume_change_diff != 0] = maxval
res_prop_voroSC$VSCmodified_volume_diff  = res_prop_voroSC$VSCmodified_volume_diff + res_prop_voroSC$VSCvolume_change_diff
res_prop_voroSC$VSCmodified_volume_ratio[res_prop_voroSC$VSCvolume_change_ratio != 1] = maxval
res_prop_voroSC$VSCmodified_volume_ratio = res_prop_voroSC$VSCmodified_volume_ratio * res_prop_voroSC$VSCvolume_change_ratio

pdb_prop_elj   = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
seqent_vorvolm_scor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( str(counter), pdb, '\n' )
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  pdb_voro = res_prop_voroSC[res_prop_voroSC$pdb==pdb,]

  x = cor.test( pdb_elj$seqent, pdb_voro$VSCvolume, method='spearman', na.action="na.omit" )
  r.seqent_vorvolSC = x$estimate
  p.seqent_vorvolSC = x$p.value
  
  x = cor.test( pdb_elj$seqent, pdb_voro$VSCmodified_volume_diff, method='spearman', na.action="na.omit" )
  r.seqent_vorvolmdSC = x$estimate
  p.seqent_vorvolmdSC = x$p.value
  
  x = cor.test( pdb_elj$seqent, pdb_voro$VSCmodified_volume_ratio, method='spearman', na.action="na.omit" )
  r.seqent_vorvolmrSC = x$estimate
  p.seqent_vorvolmrSC = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_vorvolSC,
                     r.seqent_vorvolmdSC,
                     r.seqent_vorvolmrSC
                    )

  seqent_vorvolm_scor = rbind(seqent_vorvolm_scor,srow)
}

mean(seqent_vorvolm_scor$r.seqent_vorvolSC)
mean(seqent_vorvolm_scor$r.seqent_vorvolmdSC)
mean(seqent_vorvolm_scor$r.seqent_vorvolmrSC)
median(seqent_vorvolm_scor$r.seqent_vorvolSC)
median(seqent_vorvolm_scor$r.seqent_vorvolmdSC)
median(seqent_vorvolm_scor$r.seqent_vorvolmrSC)
