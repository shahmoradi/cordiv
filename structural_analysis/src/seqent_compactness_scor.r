# Calculates the spearman correlation between Sequence Entropy (seqent) and the sphericity (compactness) of the Voronoi cell the Amino Acid.

# Amir Shahmoradi, Saturday 10:51 PM, Oct 4 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)
res_prop_voroCA      = cbind(res_prop_voroCA, VCAcompactness = sqrt(res_prop_voroCA$VCAarea^3)/res_prop_voroCA$VCAvolume)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCcompactness = sqrt(res_prop_voroSC$VSCarea^3)/res_prop_voroSC$VSCvolume)

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)
res_prop_voroAA      = cbind(res_prop_voroAA, VAAcompactness = sqrt(res_prop_voroAA$VAAarea^3)/res_prop_voroAA$VAAvolume)

pdb_prop_elj   = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
seqent_compactness_scor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( str(counter), pdb, '\n' )
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  
  pdb_voro = res_prop_voroCA[res_prop_voroCA$pdb==pdb,]
  # cat (length(pdb_voro), length(pdb_elj), '\n' )
  x = cor.test( pdb_elj$seqent, pdb_voro$VCAcompactness, method='spearman', na.action="na.omit" )
  r.seqent_compactnessCA = x$estimate
  p.seqent_compactnessCA = x$p.value
  
  pdb_voro = res_prop_voroSC[res_prop_voroSC$pdb==pdb,]
  x = cor.test( pdb_elj$seqent, pdb_voro$VSCcompactness, method='spearman', na.action="na.omit" )
  r.seqent_compactnessSC = x$estimate
  p.seqent_compactnessSC = x$p.value
  
  pdb_voro = res_prop_voroAA[res_prop_voroAA$pdb==pdb,]
  x = cor.test( pdb_elj$seqent, pdb_voro$VAAcompactness, method='spearman', na.action="na.omit" )
  r.seqent_compactnessAA = x$estimate
  p.seqent_compactnessAA = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_compactnessCA,
                     r.seqent_compactnessSC,
                     r.seqent_compactnessAA
                    )

  seqent_compactness_scor = rbind(seqent_compactness_scor,srow)
}

mean(seqent_compactness_scor$r.seqent_compactnessSC)
median(seqent_compactness_scor$r.seqent_compactnessSC)

mean(seqent_compactness_scor$r.seqent_compactnessAA)
median(seqent_compactness_scor$r.seqent_compactnessAA)

mean(seqent_compactness_scor$r.seqent_compactnessCA)
median(seqent_compactness_scor$r.seqent_compactnessCA)