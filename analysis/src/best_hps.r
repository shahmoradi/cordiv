# This code will compare the correlation strengths of the three different definitions of hydrophobicity scales with other variables, in order to find which definition is most useful in the context of globular proteins.

# Last updated by Amir Shahmoradi, Wednesday 1:58 PM, October 1 2014, Wilke Lab, ICMB, UT Austin

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

hps_scors_all_pdbs = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,c('hpskd','hpsww','hpshh')]
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ] #c('asa','rsa','hbe_mean','rss')] )
  pdb_voroAA = res_prop_voroAA[res_prop_voroAA$pdb==pdb, ]
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(seqent,ddgent)),
                    #subset(pdb_hps, select = c(hpskd,hpsww,hpshh)),
                    subset(pdb_dssp, select = c(asa,rsa,hbe_mean)),
                    subset(pdb_wcn_bf, select = -c(pdb,resnam,resnum)),
                    subset(pdb_voroAA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,VAAnvertices,VAAnedges,VAAvolume_change)),
                    subset(pdb_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change)),
                    subset(pdb_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change))
                   )
  
  pdb_hps_long = reshape(pdb_hps, ids = rownames(pdb_hps), varying = colnames(pdb_hps), v.names = 'value', timevar = 'variable', times = colnames(pdb_hps), direction = 'long')
  pdb_hps_long$variable = factor(pdb_hps_long$variable)
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  for (hps in levels(pdb_hps_long$variable))
  {
    hps_data = pdb_hps_long[pdb_hps_long$variable == hps,]
    for (variable in levels(pdb_long$variable))
    {
      variable_data = pdb_long[pdb_long$variable == variable,]
      x = cor.test( hps_data$value, variable_data$value, method='spearman', na.action="na.omit" )
      r = x$estimate
      p = x$p.value
      
      row = data.frame(pdb = pdb, hps = hps, variable = variable, value = r)
      hps_scors_all_pdbs = rbind(hps_scors_all_pdbs,row)
    }
  }
}

write.csv( hps_scors_all_pdbs, "../tables/best_HPS/hps_scors_all_pdbs.csv", row.names=F )

# Now summarize all Spearman correlations over the entire dataset

hps_scors_all_pdbs$variable = factor(hps_scors_all_pdbs$variable)

hpshh_scors_summary = data.frame()
hpskd_scors_summary = data.frame()
hpsww_scors_summary = data.frame()

for (variable in levels(hps_scors_all_pdbs$variable))
{
  cat (variable, '\n')
  temp_data = hps_scors_all_pdbs[hps_scors_all_pdbs$variable == variable,]
  if (length(temp_data$pdb) != 3*213)
  {
    cat ( 'something is fishy here!')
  }
  
  temp_data$hps = factor(temp_data$hps)
  
  # hpshh
  temp_data_hps = temp_data[temp_data$hps == 'hpshh',]
  x   = quantile(temp_data_hps$value, probs = c(0,0.25,0.5,0.75, 1.))
  row = data.frame(variable = variable,
                   mean     = mean(temp_data_hps$value),
                   median   = x[3],
                   sd       = sd(temp_data_hps$value),
                   min      = x[1],
                   q1       = x[2],
                   q3       = x[4],
                   max      = x[5]
                  )
  hpshh_scors_summary = rbind(hpshh_scors_summary,row)
  
  # hpskd
  temp_data_hps = temp_data[temp_data$hps == 'hpskd',]
  x   = quantile(temp_data_hps$value, probs = c(0,0.25,0.5,0.75, 1.))
  row = data.frame(variable = variable,
                   mean     = mean(temp_data_hps$value),
                   median   = x[3],
                   sd       = sd(temp_data_hps$value),
                   min      = x[1],
                   q1       = x[2],
                   q3       = x[4],
                   max      = x[5]
  )
  hpskd_scors_summary = rbind(hpskd_scors_summary,row)

  # hpsww
  temp_data_hps = temp_data[temp_data$hps == 'hpsww',]
  x   = quantile(temp_data_hps$value, probs = c(0,0.25,0.5,0.75, 1.))
  row = data.frame(variable = variable,
                   mean     = mean(temp_data_hps$value),
                   median   = x[3],
                   sd       = sd(temp_data_hps$value),
                   min      = x[1],
                   q1       = x[2],
                   q3       = x[4],
                   max      = x[5]
  )
  hpsww_scors_summary = rbind(hpsww_scors_summary,row)
}

row.names(hpshh_scors_summary) = c()
row.names(hpskd_scors_summary) = c()
row.names(hpsww_scors_summary) = c()
write.csv(hpshh_scors_summary, "../tables/best_HPS/hpshh_scors_summary.csv", row.names=F )
write.csv(hpskd_scors_summary, "../tables/best_HPS/hpskd_scors_summary.csv", row.names=F )
write.csv(hpsww_scors_summary, "../tables/best_HPS/hpsww_scors_summary.csv", row.names=F )

# Now compare the three HPS scors with each other:

# hpshh_hpskd:
hpshh_hpskd_comparison = data.frame(variable = hpshh_scors_summary$variable,
                                    mean_diff     = abs(hpshh_scors_summary$mean)   - abs(hpskd_scors_summary$mean),
                                    median_diff   = abs(hpshh_scors_summary$median) - abs(hpskd_scors_summary$median),
                                    sd_diff       = abs(hpshh_scors_summary$sd)     - abs(hpskd_scors_summary$sd),
                                    min_diff      = abs(hpshh_scors_summary$min)    - abs(hpskd_scors_summary$min),
                                    q1_diff       = abs(hpshh_scors_summary$q1)     - abs(hpskd_scors_summary$q1),
                                    q3_diff       = abs(hpshh_scors_summary$q3)     - abs(hpskd_scors_summary$q3),
                                    max_diff      = abs(hpshh_scors_summary$max)    - abs(hpskd_scors_summary$max)
                                    )
row.names(hpshh_hpskd_comparison) = c()
write.csv( hpshh_hpskd_comparison, "../tables/best_HPS/hpshh_hpskd_comparison.csv", row.names=F )

# hpshh_hpsww:
hpshh_hpsww_comparison = data.frame(variable = hpshh_scors_summary$variable,
                                    mean_diff     = abs(hpshh_scors_summary$mean)   - abs(hpsww_scors_summary$mean),
                                    median_diff   = abs(hpshh_scors_summary$median) - abs(hpsww_scors_summary$median),
                                    sd_diff       = abs(hpshh_scors_summary$sd)     - abs(hpsww_scors_summary$sd),
                                    min_diff      = abs(hpshh_scors_summary$min)    - abs(hpsww_scors_summary$min),
                                    q1_diff       = abs(hpshh_scors_summary$q1)     - abs(hpsww_scors_summary$q1),
                                    q3_diff       = abs(hpshh_scors_summary$q3)     - abs(hpsww_scors_summary$q3),
                                    max_diff      = abs(hpshh_scors_summary$max)    - abs(hpsww_scors_summary$max)
                                    )
row.names( hpshh_hpsww_comparison ) = c()
write.csv( hpshh_hpsww_comparison, "../tables/best_HPS/hpshh_hpsww_comparison.csv", row.names=F )

# hpskd_hpsww:
hpskd_hpsww_comparison = data.frame(variable = hpskd_scors_summary$variable,
                                    mean_diff     = abs(hpskd_scors_summary$mean)   - abs(hpsww_scors_summary$mean),
                                    median_diff   = abs(hpskd_scors_summary$median) - abs(hpsww_scors_summary$median),
                                    sd_diff       = abs(hpskd_scors_summary$sd)     - abs(hpsww_scors_summary$sd),
                                    min_diff      = abs(hpskd_scors_summary$min)    - abs(hpsww_scors_summary$min),
                                    q1_diff       = abs(hpskd_scors_summary$q1)     - abs(hpsww_scors_summary$q1),
                                    q3_diff       = abs(hpskd_scors_summary$q3)     - abs(hpsww_scors_summary$q3),
                                    max_diff      = abs(hpskd_scors_summary$max)    - abs(hpsww_scors_summary$max)
                                    )
row.names( hpskd_hpsww_comparison ) = c()
write.csv( hpskd_hpsww_comparison, "../tables/best_HPS/hpskd_hpsww_comparison.csv", row.names=F )

# Now plot and compare the correlations of the three HPS with each other.

x = -1:1

pdf( "../figures/best_HPS/hpshh_hpskd.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(hpshh_scors_summary$mean,
     hpskd_scors_summary$mean,
     xlab='Spearman r: Hydrophobicity Scale (HH)',
     ylab='Spearman r: Hydrophobicity Scale (KD)'
     )
lines( x, -x, col='red' )
graphics.off()

pdf( "../figures/best_HPS/hpshh_hpsww.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(hpshh_scors_summary$mean,
     hpsww_scors_summary$mean,
     xlab='Spearman r: Hydrophobicity Scale (HH)',
     ylab='Spearman r: Hydrophobicity Scale (WW)'
)
lines( x, -x, col='red' )
graphics.off()

pdf( "../figures/best_HPS/hpskd_hpsww.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(hpskd_scors_summary$mean,
     hpsww_scors_summary$mean,
     xlab='Spearman r: Hydrophobicity Scale (KD)',
     ylab='Spearman r: Hydrophobicity Scale (WW)'
)
lines( x, x, col='red' )
graphics.off()



# Make plots of absolute values for clarity:

x = 0:1

pdf( "../figures/best_HPS/abs(hpshh_hpskd).pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(abs(hpshh_scors_summary$mean),
     abs(hpskd_scors_summary$mean),
     xlab='Spearman r: abs(HPS HH)',
     ylab='Spearman r: abs(HPS KD)'
     )
lines( x, x, col='red' )
graphics.off()

pdf( "../figures/best_HPS/abs(hpshh_hpsww).pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(abs(hpshh_scors_summary$mean),
     abs(hpsww_scors_summary$mean),
     xlab='Spearman r: abs(HPS HH)',
     ylab='Spearman r: abs(HPS WW)'
)
lines( x, x, col='red' )
graphics.off()

pdf( "../figures/best_HPS/abs(hpskd_hpsww).pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(abs(hpskd_scors_summary$mean),
     abs(hpsww_scors_summary$mean),
     xlab='Spearman r: abs(HPS KD)',
     ylab='Spearman r: abs(HPS WW)'
)
lines( x, x, col='red' )
graphics.off()
