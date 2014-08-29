# This R function takes in the residue properties of all proteins from DSSP, and PDB files and and combines them into a single dataframe for further analysis at pdb level.

# input files:  
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out

# Last updated by Amir Shahmoradi, Thursday 3:40 PM, Aug 28 2014, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

res_prop_elj         = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb     = factor(res_prop_elj$pdb)

res_prop_hps         = read.table('../../properties/res_prop_hps.out', header=T)
res_prop_hps$pdb     = factor(res_prop_hps$pdb)

res_prop_dssp        = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp$pdb    = factor(res_prop_dssp$pdb)

res_prop_wcn_bf      = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf$pdb  = factor(res_prop_wcn_bf$pdb)

pdb_prop_long = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ] #c('asa','rsa','hbe_mean','rss')] )
  
  pdb_temp = cbind( subset(pdb_elj, select = c(seqent,ddgent)),
                    subset(pdb_hps, select = c(hpskd,hpsww,hpshh)),
                    subset(pdb_dssp, select = c(asa,rsa,hbe_mean)),
                    subset(pdb_wcn_bf, select = -c(pdb,resnam,resnum))
                  )
                  
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  for (variable1 in levels(pdb_long$variable))
  {
    #cat (variable1, '\n')
    var1 = pdb_long[pdb_long$variable == variable1,]
    
    # Claculate potentially important statistical moments of the factored variable:
      row = data.frame(pdb, variable = paste0('sum.',variable1), value = sum(var1$value))       ; pdb_prop_long = rbind(pdb_prop_long,row)
      row = data.frame(pdb, variable = paste0('mean.',variable1), value = mean(var1$value))     ; pdb_prop_long = rbind(pdb_prop_long,row)
      row = data.frame(pdb, variable = paste0('median.',variable1), value = median(var1$value)) ; pdb_prop_long = rbind(pdb_prop_long,row)
      row = data.frame(pdb, variable = paste0('sd.',variable1), value = sd(var1$value))         ; pdb_prop_long = rbind(pdb_prop_long,row)
      
    # Now calculate the Spearman correlations between pairs of variables:
      for (variable2 in levels(pdb_long$variable))
      {
        if ( variable1 != variable2)
        {
          var2 = pdb_long[pdb_long$variable == variable2,]
          x = cor.test( var1$value, var2$value, method='spearman', na.action="na.omit" )
          r = x$estimate
          p = x$p.value
          
          row = data.frame(pdb, variable = paste0('r.',variable1,'.',variable2), value = r)
          pdb_prop_long = rbind(pdb_prop_long,row)
        }
      }
  }

}


  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$asa, method='spearman', na.action="na.omit" )
  r.seqent_asa = x$estimate
  p.seqent_asa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$rsa, method='spearman', na.action="na.omit" )
  r.seqent_rsa = x$estimate
  p.seqent_rsa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$wcn_ca, method='spearman', na.action="na.omit" )
  r.seqent_wcnca = x$estimate
  p.seqent_wcnca = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_alignments, pdb_dssp$hbe_mean, method='spearman', na.action="na.omit" )
  r.seqent_hbe = x$estimate
  p.seqent_hbe = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$asa, method='spearman', na.action="na.omit" )
  r.ddgent_asa = x$estimate
  p.ddgent_asa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$rsa, method='spearman', na.action="na.omit" )
  r.ddgent_rsa = x$estimate
  p.ddgent_rsa = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$wcn_ca, method='spearman', na.action="na.omit" )
  r.ddgent_wcnca = x$estimate
  p.ddgent_wcnca = x$p.value
  
  x = cor.test( pdb_elj$entropy_from_ddGs, pdb_dssp$hbe_mean, method='spearman', na.action="na.omit" )
  r.ddgent_hbe = x$estimate
  p.ddgent_hbe = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$rsa, method='spearman', na.action="na.omit" )
  r.asa_rsa = x$estimate
  p.asa_rsa = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$wcn_ca, method='spearman', na.action="na.omit" )
  r.asa_wcnca = x$estimate
  p.asa_wcnca = x$p.value
  
  x = cor.test( pdb_dssp$asa, pdb_dssp$hbe_mean, method='spearman', na.action="na.omit" )
  r.asa_hbe = x$estimate
  p.asa_hbe = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_dssp$wcn_ca, method='spearman', na.action="na.omit" )
  r.rsa_wcnca = x$estimate
  p.rsa_wcnca = x$p.value
  
  x = cor.test( pdb_dssp$rsa, pdb_dssp$hbe_mean, method='spearman', na.action="na.omit" )
  r.rsa_hbe = x$estimate
  p.rsa_hbe = x$p.value
  
  x = cor.test( pdb_dssp$wcn_ca, pdb_dssp$hbe_mean, method='spearman', na.action="na.omit" )
  r.wcnca_hbe = x$estimate
  p.wcnca_hbe = x$p.value
  
  srow = data.frame( pdb=pdb,
                     r.seqent_ddgent = r.seqent_ddgent,
                     r.seqent_asa    = r.seqent_asa,
                     r.seqent_rsa    = r.seqent_rsa,    
                     r.seqent_wcnca  = r.seqent_wcnca,
                     r.seqent_hbe    = r.seqent_hbe,    
                     r.ddgent_asa    = r.ddgent_asa,    
                     r.ddgent_rsa    = r.ddgent_rsa,    
                     r.ddgent_wcnca  = r.ddgent_wcnca,
                     r.ddgent_hbe    = r.ddgent_hbe,    
                     r.asa_rsa       = r.asa_rsa,       
                     r.asa_wcnca     = r.asa_wcnca,
                     r.asa_hbe       = r.asa_hbe,       
                     r.rsa_wcnca     = r.rsa_wcnca,
                     r.rsa_hbe       = r.rsa_hbe,       
                     r.wcnca_hbe     = r.wcnca_hbe  )

  prow = data.frame( pdb=pdb,
                     p.seqent_ddgent = p.seqent_ddgent,
                     p.seqent_asa    = p.seqent_asa,    
                     p.seqent_rsa    = p.seqent_rsa,    
                     p.seqent_wcnca  = p.seqent_wcnca,  
                     p.seqent_hbe    = p.seqent_hbe,    
                     p.ddgent_asa    = p.ddgent_asa,    
                     p.ddgent_rsa    = p.ddgent_rsa,    
                     p.ddgent_wcnca  = p.ddgent_wcnca,  
                     p.ddgent_hbe    = p.ddgent_hbe,    
                     p.asa_rsa       = p.asa_rsa,       
                     p.asa_wcnca     = p.asa_wcnca,
                     p.asa_hbe       = p.asa_hbe,       
                     p.rsa_wcnca     = p.rsa_wcnca,
                     p.rsa_hbe       = p.rsa_hbe,       
                     p.wcnca_hbe     = p.wcnca_hbe  )
  
  erow = data.frame( pdb=pdb,
                     sum.seqent    = sum(pdb_elj$entropy_from_alignments),
                     mean.seqent   = mean(pdb_elj$entropy_from_alignments),
                     median.seqent = median(pdb_elj$entropy_from_alignments),
                     sd.seqent     = sd(pdb_elj$entropy_from_alignments),
                     sum.ddgent    = sum(pdb_elj$entropy_from_ddGs),
                     mean.ddgent   = mean(pdb_elj$entropy_from_ddGs),
                     median.ddgent = median(pdb_elj$entropy_from_ddGs),
                     sd.ddgent     = sd(pdb_elj$entropy_from_ddGs)   )
  
  pdb_prop_scor  = rbind(pdb_prop_scor,srow)
  pdb_prop_scorp = rbind(pdb_prop_scorp,prow)
  pdb_prop_elj   = rbind(pdb_prop_elj,erow)

}

row.names(pdb_prop_scor)  = c()
row.names(pdb_prop_scorp) = c()
row.names(pdb_prop_elj) = c()

write.csv( pdb_prop_scor, "../tables/pdb_prop_scor.csv", row.names=F )
write.csv( pdb_prop_scorp, "../tables/pdb_prop_scorp.csv", row.names=F )

all_pdb_prop = cbind(pdb_prop_dssp,subset(pdb_prop_scor, select = -c(pdb)),subset(pdb_prop_elj, select = -c(pdb)))

all_pdb_prop_subset = subset(all_pdb_prop, select = -c(name,nssb))
cormat = cor(all_pdb_prop_subset, method = 'spearman')
#corrplot.mixed(cormat)
corrplot(cormat, method='circle')

# Now do some statistics on the calculated correlations. First melt the correlation data.frame to make a long format data set.
pdb_prop_scor           = melt(pdb_prop_scor)
pdb_prop_scor$variable  = factor(pdb_prop_scor$variable)

scor_stat = data.frame()

for (variable in levels(pdb_prop_scor$variable))
{
  temp = pdb_prop_scor[pdb_prop_scor$variable==variable,]
  mean.scor   = mean(temp$value)
  median.scor = median(temp$value)
  sd.scor     = sd(temp$value)
  
  row = data.frame(variable = variable,
                   mean     = mean.scor,
                   median   = median.scor,
                   sd       = sd.scor
                   )
  scor_stat = rbind(scor_stat,row)
}

row.names(scor_stat) = c()
write.csv( scor_stat, "../tables/scor_stat.csv", row.names=F )
