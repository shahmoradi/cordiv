# This R script reads in Residue data WCN and Bfactor (BF) for different atoms of the residue and tests which atom best represnt WCN in predicting sequence evariability.

# Amir Shahmoradi, Thusday 4:20 PM, Aug 14 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_dssp      = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp$pdb  = factor(res_prop_dssp$pdb)

res_prop_wcn_bf = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf$pdb = factor(res_prop_wcn_bf$pdb)

pdb_prop_seqent_wcn_scor  = data.frame()

counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  pdb_wcn  = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb,]
  pdb_dssp = res_prop_dssp[res_prop_dssp$pdb==pdb,]
  
  # correlations with seqent with WCN

  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnSC, method='spearman', na.action="na.omit" )
  r.seqent_wcnSC = x$estimate
  p.seqent_wcnSC = x$p.value

  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnAA, method='spearman', na.action="na.omit" )
  r.seqent_wcnAA = x$estimate
  p.seqent_wcnAA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnN, method='spearman', na.action="na.omit" )
  r.seqent_wcnN = x$estimate
  p.seqent_wcnN = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnCA, method='spearman', na.action="na.omit" )
  r.seqent_wcnCA = x$estimate
  p.seqent_wcnCA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnC, method='spearman', na.action="na.omit" )
  r.seqent_wcnC = x$estimate
  p.seqent_wcnC = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnO, method='spearman', na.action="na.omit" )
  r.seqent_wcnO = x$estimate
  p.seqent_wcnO = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$wcnCB, method='spearman', na.action="na.omit" )
  r.seqent_wcnCB = x$estimate
  p.seqent_wcnCB = x$p.value
  
  # Now correlations of seqent with Bfactors
    
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfSC, method='spearman', na.action="na.omit" )
  r.seqent_bfSC = x$estimate
  p.seqent_bfSC = x$p.value

  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfAA, method='spearman', na.action="na.omit" )
  r.seqent_bfAA = x$estimate
  p.seqent_bfAA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfN, method='spearman', na.action="na.omit" )
  r.seqent_bfN = x$estimate
  p.seqent_bfN = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfCA, method='spearman', na.action="na.omit" )
  r.seqent_bfCA = x$estimate
  p.seqent_bfCA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfC, method='spearman', na.action="na.omit" )
  r.seqent_bfC = x$estimate
  p.seqent_bfC = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfO, method='spearman', na.action="na.omit" )
  r.seqent_bfO = x$estimate
  p.seqent_bfO = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_wcn$bfCB, method='spearman', na.action="na.omit" )
  r.seqent_bfCB = x$estimate
  p.seqent_bfCB = x$p.value
  
  
  # correlations with ddgent
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnSC, method='spearman', na.action="na.omit" )
  r.ddgent_wcnSC = x$estimate
  p.ddgent_wcnSC = x$p.value

  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnAA, method='spearman', na.action="na.omit" )
  r.ddgent_wcnAA = x$estimate
  p.ddgent_wcnAA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnN, method='spearman', na.action="na.omit" )
  r.ddgent_wcnN = x$estimate
  p.ddgent_wcnN = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnCA, method='spearman', na.action="na.omit" )
  r.ddgent_wcnCA = x$estimate
  p.ddgent_wcnCA = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnC, method='spearman', na.action="na.omit" )
  r.ddgent_wcnC = x$estimate
  p.ddgent_wcnC = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnO, method='spearman', na.action="na.omit" )
  r.ddgent_wcnO = x$estimate
  p.ddgent_wcnO = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_wcn$wcnCB, method='spearman', na.action="na.omit" )
  r.ddgent_wcnCB = x$estimate
  p.ddgent_wcnCB = x$p.value
  
  # correlations with RSA
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnSC, method='spearman', na.action="na.omit" )
  r.rsa_wcnSC = x$estimate
  p.rsa_wcnSC = x$p.value

  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnAA, method='spearman', na.action="na.omit" )
  r.rsa_wcnAA = x$estimate
  p.rsa_wcnAA = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnN, method='spearman', na.action="na.omit" )
  r.rsa_wcnN = x$estimate
  p.rsa_wcnN = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnCA, method='spearman', na.action="na.omit" )
  r.rsa_wcnCA = x$estimate
  p.rsa_wcnCA = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnC, method='spearman', na.action="na.omit" )
  r.rsa_wcnC = x$estimate
  p.rsa_wcnC = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnO, method='spearman', na.action="na.omit" )
  r.rsa_wcnO = x$estimate
  p.rsa_wcnO = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_wcn$wcnCB, method='spearman', na.action="na.omit" )
  r.rsa_wcnCB = x$estimate
  p.rsa_wcnCB = x$p.value
  
  
  srow = data.frame( pdb=pdb,

                     r.seqent_wcnSC = r.seqent_wcnSC,
                     r.seqent_wcnAA = r.seqent_wcnAA,
                     r.seqent_wcnN  = r.seqent_wcnN,
                     r.seqent_wcnCA = r.seqent_wcnCA,
                     r.seqent_wcnC  = r.seqent_wcnC,
                     r.seqent_wcnO  = r.seqent_wcnO,
                     r.seqent_wcnCB = r.seqent_wcnCB,

                     r.seqent_bfSC = r.seqent_bfSC,
                     r.seqent_bfAA = r.seqent_bfAA,
                     r.seqent_bfN  = r.seqent_bfN,
                     r.seqent_bfCA = r.seqent_bfCA,
                     r.seqent_bfC  = r.seqent_bfC,
                     r.seqent_bfO  = r.seqent_bfO,
                     r.seqent_bfCB = r.seqent_bfCB,

                     r.ddgent_wcnSC = r.ddgent_wcnSC,
                     r.ddgent_wcnAA = r.ddgent_wcnAA,
                     r.ddgent_wcnN  = r.ddgent_wcnN,
                     r.ddgent_wcnCA = r.ddgent_wcnCA,
                     r.ddgent_wcnC  = r.ddgent_wcnC,
                     r.ddgent_wcnO  = r.ddgent_wcnO,
                     r.ddgent_wcnCB = r.ddgent_wcnCB,

                     r.rsa_wcnSC = r.rsa_wcnSC,
                     r.rsa_wcnAA = r.rsa_wcnAA,
                     r.rsa_wcnN  = r.rsa_wcnN,
                     r.rsa_wcnCA = r.rsa_wcnCA,
                     r.rsa_wcnC  = r.rsa_wcnC,
                     r.rsa_wcnO  = r.rsa_wcnO,
                     r.rsa_wcnCB = r.rsa_wcnCB
                     
                     )
                     
  pdb_prop_seqent_wcn_scor = rbind(pdb_prop_seqent_wcn_scor,srow)
  
  cat( str(counter), pdb )

}

row.names(pdb_prop_seqent_wcn_scor)  = c()

write.csv( pdb_prop_seqent_wcn_scor, "../tables/pdb_prop_seqent_wcn_scor.csv", row.names=F )
