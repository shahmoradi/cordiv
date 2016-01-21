# This sctipt compares the calculated rates from the two evolutioray models JTT and JC

JC_JTT_spcor = data.frame()
counter = 0

# Rate4Site site-specific evolutionary rates
res_prop_jec         = read.csv('../../jec_pdb_r4s.csv', header=T)
res_prop_jec         = res_prop_jec[!(res_prop_jec$pdb %in% excluded_pdbs),]
res_prop_jec$pdb     = factor(res_prop_jec$pdb)

for (pdb in levels(res_prop_jec$pdb))
{
  counter = counter + 1
  print( paste (counter, ":" , pdb) )
  data = res_prop_jec[res_prop_jec$pdb == pdb,]
  x = cor.test( data$r4s_JTT, data$r4s_JC, method='spearman', na.action="na.omit" )
  r.JC_JTT = x$estimate
  p.JC_JTT = x$p.value
  row = data.frame( pdb=pdb,
                    r.JC_JTT = r.JC_JTT
                  , p.JC_JTT = p.JC_JTT
  )
  JC_JTT_spcor = rbind(JC_JTT_spcor,row)                      
}

row.names(JC_JTT_spcor) = c()
write.csv( JC_JTT_spcor, "../tables/JC_JTT_spcor.csv", row.names=F )

mean(JC_JTT_spcor$r.JC_JTT)
sd(JC_JTT_spcor$r.JC_JTT)
