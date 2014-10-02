# I am writing this code to find out how strong different residue properties correlate with each other, when averaged over all PDB dataset.
# The output will be a table that contains, on each row, the mean, median and STDEV of the Spearman correlations of pairs of residue properties over the entire PDB dataset.

# Input data:
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out
#               ../../properties/res_prop_voronoiAA.out
#               ../../properties/res_prop_voronoiCA.out
#               ../../properties/res_prop_voronoiSC.out

# Last updated by Amir Shahmoradi, Tuesday 7:41 PM, Sep 30 2014, Wilke Lab, ICMB, UT Austin

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

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)
res_prop_voroAA      = cbind(res_prop_voroAA, VAAmodified_volume = res_prop_voroAA$VAAvolume)
maxval = max(res_prop_voroAA$VAAvolume)
res_prop_voroAA$VAAmodified_volume[res_prop_voroAA$VAAvolume_change != 0] = maxval
res_prop_voroAA$VAAmodified_volume = res_prop_voroAA$VAAmodified_volume + res_prop_voroAA$VAAvolume_change

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)
res_prop_voroCA      = cbind(res_prop_voroCA, VCAmodified_volume = res_prop_voroCA$VCAvolume)
maxval = max(res_prop_voroCA$VCAvolume)
res_prop_voroCA$VCAmodified_volume[res_prop_voroCA$VCAvolume_change != 0] = maxval
res_prop_voroCA$VCAmodified_volume = res_prop_voroCA$VCAmodified_volume + res_prop_voroCA$VCAvolume_change

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCmodified_volume = res_prop_voroSC$VSCvolume)
maxval = max(res_prop_voroSC$VSCvolume)
res_prop_voroSC$VSCmodified_volume[res_prop_voroSC$VSCvolume_change != 0] = maxval
res_prop_voroSC$VSCmodified_volume = res_prop_voroSC$VSCmodified_volume + res_prop_voroSC$VSCvolume_change

all_scors_all_pdbs = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ] #c('asa','rsa','hbe_mean','rss')] )
  pdb_voroAA = res_prop_voroAA[res_prop_voroAA$pdb==pdb, ]
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(seqent,ddgent)),
                    subset(pdb_hps, select = c(hpskd,hpsww,hpshh)),
                    subset(pdb_dssp, select = c(asa,rsa,hbe_mean)),
                    subset(pdb_wcn_bf, select = -c(pdb,resnam,resnum)),
                    subset(pdb_voroAA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,VAAnvertices,VAAnedges,VAAvolume_change)),
                    subset(pdb_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change)),
                    subset(pdb_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change))
                  )
                  
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  counter1 = 0
  
  for (variable1 in levels(pdb_long$variable))
  {
    counter1 = counter1 + 1
    #cat (variable1, '\n')
    var1 = pdb_long[pdb_long$variable == variable1,]
    
    # Now calculate the Spearman correlations between pairs of variables:
    counter2 = 0
    for (variable2 in levels(pdb_long$variable))
    {
      counter2 = counter2 + 1
      if ( variable1 != variable2 & counter1 < counter2)
      { 
        # cat ( variable1,variable2,str(length(variable1)),str(length(variable2)) )
        var2 = pdb_long[pdb_long$variable == variable2,]
        x = cor.test( var1$value, var2$value, method='spearman', na.action="na.omit" )
        r = x$estimate
        p = x$p.value
        
        row = data.frame(pdb = pdb, variable = paste0('r.',variable1,'.',variable2), value = r)
        all_scors_all_pdbs = rbind(all_scors_all_pdbs,row)
      }
    }
  }

}

write.csv( all_scors_all_pdbs, "../tables/all_scors_all_pdbs.csv", row.names=F )

# Now summarize all Spearman correlations over the entire dataset

all_scors_all_pdbs$variable = factor(all_scors_all_pdbs$variable)

all_scors_summary = data.frame()

for (variable in levels(all_scors_all_pdbs$variable))
{
  cat (variable, '\n')
  temp_data = all_scors_all_pdbs[all_scors_all_pdbs$variable == variable,]
  if (length(temp_data$pdb) != 213)
  {
    cat ( 'something is fishy here!')
  }

  x = quantile(temp_data$value, probs = c(0,0.25,0.5,0.75, 1.))
  
  row = data.frame(variable = variable,
                   mean     = mean(temp_data$value),
                   median   = x[3],
                   sd       = sd(temp_data$value),
                   min      = x[1],
                   q1       = x[2],
                   q3       = x[4],
                   max      = x[5]
                   )
  all_scors_summary = rbind(all_scors_summary,row)
}

row.names(all_scors_summary) = c()
write.csv(all_scors_summary, "../tables/all_scors_summary.csv", row.names=F )

