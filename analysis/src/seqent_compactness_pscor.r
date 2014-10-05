# Calculates the PARTIAL spearman correlation between Sequence Entropy (seqent) and the sphericity (sphericity) of the Voronoi cell the Amino Acid, controlling for Voronoi Volume

# Amir Shahmoradi, Saturday 10:51 PM, Oct 4 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)
res_prop_voroCA      = cbind(res_prop_voroCA, VCAsphericity = 4.8359758620494089221509005399179*(res_prop_voroCA$VCAvolume^(2./3.))/res_prop_voroCA$VCAarea)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC$VSCvolume^(2./3.))/res_prop_voroSC$VSCarea)

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)
res_prop_voroAA      = cbind(res_prop_voroAA, VAAsphericity = 4.8359758620494089221509005399179*(res_prop_voroAA$VAAvolume^(2./3.))/res_prop_voroAA$VAAarea)

pdb_prop_elj   = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
seqent_sphericity_pscor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( str(counter), pdb, '\n' )
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb,]
  # cat (length(pdb_voroCA), length(pdb_elj), '\n' )
  x = pcor.test( pdb_elj$seqent, pdb_voroCA$VCAsphericity, pdb_voroCA$VCAvolume, method='spearman' )
  r.seqent_sphericityCA = x$estimate
  p.seqent_sphericityCA = x$p.value
  
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb,]
  x = pcor.test( pdb_elj$seqent, pdb_voroSC$VSCsphericity, pdb_voroSC$VSCvolume, method='spearman' )
  r.seqent_sphericitySC = x$estimate
  p.seqent_sphericitySC = x$p.value
  
  pdb_voroAA = res_prop_voroAA[res_prop_voroAA$pdb==pdb,]
  x = pcor.test( pdb_elj$seqent, pdb_voroAA$VAAsphericity, pdb_voroAA$VAAvolume, method='spearman' )
  r.seqent_sphericityAA = x$estimate
  p.seqent_sphericityAA = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_sphericityCA,
                     r.seqent_sphericitySC,
                     r.seqent_sphericityAA
                    )

  seqent_sphericity_pscor = rbind(seqent_sphericity_pscor,srow)
}

mean(seqent_sphericity_pscor$r.seqent_sphericitySC)
median(seqent_sphericity_pscor$r.seqent_sphericitySC)

mean(seqent_sphericity_pscor$r.seqent_sphericityAA)
median(seqent_sphericity_pscor$r.seqent_sphericityAA)

mean(seqent_sphericity_pscor$r.seqent_sphericityCA)
median(seqent_sphericity_pscor$r.seqent_sphericityCA)