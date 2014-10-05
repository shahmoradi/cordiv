# Calculates the spearman correlation between Sequence Entropy (seqent) and the Voronoi surface areas of the Amino Acid.

# Amir Shahmoradi, Friday 1:25 PM, Sep 19 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

pdb_prop_elj   = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
seqent_vorarea_scor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( str(counter), pdb, '\n' )
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  
  pdb_voro = res_prop_voroCA[res_prop_voroCA$pdb==pdb,]
  # cat (length(pdb_voro), length(pdb_elj), '\n' )
  x = cor.test( pdb_elj$seqent, pdb_voro$VCAarea, method='spearman', na.action="na.omit" )
  r.seqent_vorareaCA = x$estimate
  p.seqent_vorareaCA = x$p.value
  
  pdb_voro = res_prop_voroSC[res_prop_voroSC$pdb==pdb,]
  x = cor.test( pdb_elj$seqent, pdb_voro$VSCarea, method='spearman', na.action="na.omit" )
  r.seqent_vorareaSC = x$estimate
  p.seqent_vorareaSC = x$p.value
  
  pdb_voro = res_prop_voroAA[res_prop_voroAA$pdb==pdb,]
  x = cor.test( pdb_elj$seqent, pdb_voro$VAAarea, method='spearman', na.action="na.omit" )
  r.seqent_vorareaAA = x$estimate
  p.seqent_vorareaAA = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_vorareaCA,
                     r.seqent_vorareaSC,
                     r.seqent_vorareaAA
                    )

  seqent_vorarea_scor = rbind(seqent_vorarea_scor,srow)
}

mean(seqent_vorarea_scor$r.seqent_vorareaSC)
median(seqent_vorarea_scor$r.seqent_vorareaSC)

mean(seqent_vorarea_scor$r.seqent_vorareaAA)
median(seqent_vorarea_scor$r.seqent_vorareaAA)

mean(seqent_vorarea_scor$r.seqent_vorareaCA)
median(seqent_vorarea_scor$r.seqent_vorareaCA)