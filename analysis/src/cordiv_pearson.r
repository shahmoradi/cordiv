# This R function was written for the purpose of finding out the factors affecting the strength of the correlation between structural properties and sequence evolution.
# There three primary input files containing data:  elj_pdb_entropies.in  +  res_prop_dssp.out + pdb_properties.out

# Amir Shahmoradi, Monday 3:20 PM, Aug 4 2014, Wilke Lab, ICMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

res_prop_elj       = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb   = factor(res_prop_elj$pdb)

res_prop_dssp      = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp$pdb  = factor(res_prop_dssp$pdb)

pdb_prop_dssp = read.table('../../properties/pdb_prop_dssp.out', header=T)

pdb_prop_pcor = data.frame()
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  pdb_elj  = res_prop_elj[res_prop_elj$pdb==pdb,]
  pdb_dssp = res_prop_dssp[res_prop_dssp$pdb==pdb,]
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_elj$entropy_from_ddGs, method='pearson', na.action="na.omit" )
  r.seqent_ddgent = x$estimate
  p.seqent_ddgent = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$asa, method='pearson', na.action="na.omit" )
  r.seqent_asa = x$estimate
  p.seqent_asa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$rsa, method='pearson', na.action="na.omit" )
  r.seqent_rsa = x$estimate
  p.seqent_rsa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$wcn_ca, method='pearson', na.action="na.omit" )
  r.seqent_wcnca = x$estimate
  p.seqent_wcnca = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$hbe_mean, method='pearson', na.action="na.omit" )
  r.seqent_hbe = x$estimate
  p.seqent_hbe = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$asa, method='pearson', na.action="na.omit" )
  r.ddgent_asa = x$estimate
  p.ddgent_asa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$rsa, method='pearson', na.action="na.omit" )
  r.ddgent_rsa = x$estimate
  p.ddgent_rsa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$wcn_ca, method='pearson', na.action="na.omit" )
  r.ddgent_wcnca = x$estimate
  p.ddgent_wcnca = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$hbe_mean, method='pearson', na.action="na.omit" )
  r.ddgent_hbe = x$estimate
  p.ddgent_hbe = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$rsa, method='pearson', na.action="na.omit" )
  r.asa_rsa = x$estimate
  p.asa_rsa = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$wcn_ca, method='pearson', na.action="na.omit" )
  r.asa_wcnca = x$estimate
  p.asa_wcnca = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$hbe_mean, method='pearson', na.action="na.omit" )
  r.asa_hbe = x$estimate
  p.asa_hbe = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_dssp$wcn_ca, method='pearson', na.action="na.omit" )
  r.rsa_wcnca = x$estimate
  p.rsa_wcnca = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_dssp$hbe_mean, method='pearson', na.action="na.omit" )
  r.rsa_hbe = x$estimate
  p.rsa_hbe = x$p.value
  
  x = cor.test( pdb_dssp$wcn_ca, pdb_dssp$hbe_mean, method='pearson', na.action="na.omit" )
  r.wcnca_hbe = x$estimate
  p.wcnca_hbe = x$p.value
  
  row = data.frame( pdb=pdb,
                    r.seqent_ddgent = r.seqent_ddgent,
                    p.seqent_ddgent = p.seqent_ddgent,
                    r.seqent_asa    = r.seqent_asa,
                    p.seqent_asa    = p.seqent_asa,    
                    r.seqent_rsa    = r.seqent_rsa,    
                    p.seqent_rsa    = p.seqent_rsa,    
                    r.seqent_wcnca  = r.seqent_wcnca,
                    p.seqent_wcnca  = p.seqent_wcnca,  
                    r.seqent_hbe    = r.seqent_hbe,    
                    p.seqent_hbe    = p.seqent_hbe,    
                    r.ddgent_asa    = r.ddgent_asa,    
                    p.ddgent_asa    = p.ddgent_asa,    
                    r.ddgent_rsa    = r.ddgent_rsa,    
                    p.ddgent_rsa    = p.ddgent_rsa,    
                    r.ddgent_wcnca  = r.ddgent_wcnca,
                    p.ddgent_wcnca  = p.ddgent_wcnca,  
                    r.ddgent_hbe    = r.ddgent_hbe,    
                    p.ddgent_hbe    = p.ddgent_hbe,    
                    r.asa_rsa       = r.asa_rsa,       
                    p.asa_rsa       = p.asa_rsa,       
                    r.asa_wcnca     = r.asa_wcnca,
                    p.asa_wcnca     = p.asa_wcnca,
                    r.asa_hbe       = r.asa_hbe,       
                    p.asa_hbe       = p.asa_hbe,       
                    r.rsa_wcnca     = r.rsa_wcnca,
                    p.rsa_wcnca     = p.rsa_wcnca,
                    r.rsa_hbe       = r.rsa_hbe,       
                    p.rsa_hbe       = p.rsa_hbe,       
                    r.wcnca_hbe     = r.wcnca_hbe,     
                    p.wcnca_hbe     = p.wcnca_hbe  )
                    
  
  pdb_prop_pcor = rbind(pdb_prop_pcor,row)
  print( pdb )
  
}

row.names(pdb_prop_pcor) = c()
write.csv( pdb_prop_pcor, "../tables/pdb_prop_pcor.csv", row.names=F )