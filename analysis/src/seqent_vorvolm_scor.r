# Calculates the spearman correlation between Sequence Entropy (seqent) and the MODIFIED Voronoi volumes of the Amino Acid.

# Amir Shahmoradi, Friday 5:58 PM, Sep 19 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

#res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
#res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

#res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
#res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
res_prop_voroSC      = cbind(res_prop_voroSC, modified_volume = res_prop_voroSC$volume)
maxval = max(res_prop_voroSC$volume)
res_prop_voroSC$modified_volume[res_prop_voroSC$volume_change != 0] = maxval
res_prop_voroSC$modified_volume = res_prop_voroSC$modified_volume + res_prop_voroSC$volume_change
res_prop_voroSC      = cbind(res_prop_voroSC, modified_volume_normed_sizeSC = res_prop_voroSC$modified_volume/res_prop_voroSC$sizeSC)

pdb_prop_elj   = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
seqent_vorvolm_scor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( str(counter), pdb, '\n' )
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  pdb_voro = res_prop_voroSC[res_prop_voroSC$pdb==pdb,]

  x = cor.test( pdb_elj$seqent, pdb_voro$modified_volume, method='spearman', na.action="na.omit" )
  r.seqent_vorvolmSC = x$estimate
  p.seqent_vorvolmSC = x$p.value
  
  x = cor.test( pdb_elj$seqent, pdb_voro$modified_volume_normed_sizeSC, method='spearman', na.action="na.omit" )
  r.seqent_vorvolmnSC = x$estimate
  p.seqent_vorvolmnSC = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_vorvolmSC,
                     r.seqent_vorvolmnSC
                    )

  seqent_vorvolm_scor = rbind(seqent_vorvolm_scor,srow)
}

mean(seqent_vorvolm_scor$r.seqent_vorvolmSC)
mean(seqent_vorvolm_scor$r.seqent_vorvolmnSC)
median(seqent_vorvolm_scor$r.seqent_vorvolmSC)
median(seqent_vorvolm_scor$r.seqent_vorvolmnSC)
