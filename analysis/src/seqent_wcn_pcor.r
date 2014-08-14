# This R script reads in Residue data WCN and Bfactor (BF) for different atoms of the residue and tests which atom best represnt WCN in predicting sequence evariability.
# The only difference of this code with seqent_wcn_scor.r is that here the Pearosn correlations are calculated and output.

# Amir Shahmoradi, Thusday 4:40 PM, Aug 14 2014, Wilke Lab, ICMB, UT Austin

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_wcn_bf = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf$pdb = factor(res_prop_wcn_bf$pdb)

pdb_prop_seqent_wcn_pcor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  pdb_wcn  = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb,]
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnSC, method='pearson', na.action="na.omit" )
  r.seqent_wcnSC = x$estimate
  p.seqent_wcnSC = x$p.value

  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnAA, method='pearson', na.action="na.omit" )
  r.seqent_wcnAA = x$estimate
  p.seqent_wcnAA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnN, method='pearson', na.action="na.omit" )
  r.seqent_wcnN = x$estimate
  p.seqent_wcnN = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnCA, method='pearson', na.action="na.omit" )
  r.seqent_wcnCA = x$estimate
  p.seqent_wcnCA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnC, method='pearson', na.action="na.omit" )
  r.seqent_wcnC = x$estimate
  p.seqent_wcnC = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnO, method='pearson', na.action="na.omit" )
  r.seqent_wcnO = x$estimate
  p.seqent_wcnO = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnCB, method='pearson', na.action="na.omit" )
  r.seqent_wcnCB = x$estimate
  p.seqent_wcnCB = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_wcnSC = r.seqent_wcnSC,
                     r.seqent_wcnAA = r.seqent_wcnAA,
                     r.seqent_wcnN  = r.seqent_wcnN,
                     r.seqent_wcnCA = r.seqent_wcnCA,
                     r.seqent_wcnC  = r.seqent_wcnC,
                     r.seqent_wcnO  = r.seqent_wcnO,
                     r.seqent_wcnCB = r.seqent_wcnCB)
                     
  pdb_prop_seqent_wcn_pcor = rbind(pdb_prop_seqent_wcn_pcor,srow)
  
  cat( str(counter), pdb )

}

row.names(pdb_prop_seqent_wcn_pcor)  = c()

write.csv( pdb_prop_seqent_wcn_pcor, "../tables/pdb_prop_seqent_wcn_pcor.csv", row.names=F )