# Amir Shahmoradi, Wednesday 3:23 PM, Sep 24 2014, Wilke Lab, ICMB, UT Austin

# install.packages('zoo')
library('zoo')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#volume_smooth = filter(res_prop_all_ordered$volume, ma1000)
#plot(res_prop_all_ordered$wcnSC,volume_smooth, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )
#plot(smoothed$wcnSC_mean,smoothed$vol_mean, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )

# The following is short version of the most useful residue properties that I have found so far.
res_prop_concise = data.frame(seqent        = res_prop_elj$seqent,
                              ddgent        = res_prop_elj$ddgent,
                              rsa           = res_prop_dssp$rsa,
                              hbe           = res_prop_dssp$hbe_mean,
                              hpshh         = res_prop_hps$hpshh,
                              wcnSC         = res_prop_wcn_bf$wcnSC,
                              bfSC          = res_prop_wcn_bf$bfSC,
                              #resvol        = res_prop_voroSC$resvol,
                              #vnvertices    = res_prop_voroSC$VSCnvertices,
                              #vnedges       = res_prop_voroSC$VSCnedges,
                              vnfaces       = res_prop_voroSC$VSCnfaces,
                              vedge         = res_prop_voroSC$VSCedge_length_total,
                              varea         = res_prop_voroSC$VSCarea,
                              vvolume       = res_prop_voroSC$VSCvolume,
                              veccentricity = res_prop_voroSC$VSCeccentricity,
                              vsphericity   = res_prop_voroSC$VSCsphericity
                              )

cormat = cor(res_prop_concise, method='spearman')
write.csv( cormat, "../tables/res_prop_cormat.csv", row.names=T )

# Now since the three Voronoi quantities vnvortices, vnedges and vnfaces happen to be exactly the same as seen in the cormat calculated above, I am going to remove them from the data set, in addition to resvol which is not as informative.
# res_prop_concise = subset(res_prop_concise, select = -c(vnvertices,vnedges,resvol))
res_prop_concise_closed = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff == 0,]
res_prop_concise_open = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff != 0,]
# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Sequence Entropy' , 'ddG Entropy' , 'Relative Solvent Accessibility' , 'Hydrogen Bond Energy' ,
                  'Hydrophobicity Scale (HH)' , 'Side-Chain Contact Number' ,
                  'Average Side-Chain B-factor' , 'Voronoi Cell Faces' , 'Voronoi Cell Edge length' , 'Voronoi Cell Surface Area' ,
                  'Voronoi Cell Volume' , 'Voronoi Cell Eccentricity' , 'Voronoi cell Sphericity')

varnames_short = colnames(res_prop_concise)

for (i in 1:length(varnames_short))
{
  cat(i, varnames_short[i], '\n')
  
  res_prop_ordered = res_prop_concise[with(res_prop_concise, order(res_prop_concise[[varnames_short[[i]][1]]])),]
  test = rollapply(res_prop_ordered, width = 2000, FUN = mean)
  test = data.frame(test)
  for (j in 1:length(varnames_short))
  {
    if (j != i)
    {
      filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_all_cells.pdf')
      pdf( filename, width=4.5, height=4, useDingbats=FALSE )
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      plot(test[[varnames_short[[i]][1]]],
           test[[varnames_short[[j]][1]]],
           xlab = varnames_long[i],
           ylab = varnames_long[j],
           type='l')
      graphics.off()
    }
  }
  
  res_prop_ordered_closed = res_prop_concise_closed[with(res_prop_concise_closed, order(res_prop_concise_closed[[varnames_short[[i]][1]]])),]
  test = rollapply(res_prop_ordered_closed, width = 2000, FUN = mean)
  test = data.frame(test)
  for (j in 1:length(varnames_short))
  {
    if (j != i)
    {
      filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_closed_cells.pdf')
      pdf( filename, width=4.5, height=4, useDingbats=FALSE )
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      plot(test[[varnames_short[[i]][1]]],
           test[[varnames_short[[j]][1]]],
           xlab = varnames_long[i],
           ylab = varnames_long[j],
           type='l')
      graphics.off()
    }
  }
  
  res_prop_ordered_open = res_prop_concise_open[with(res_prop_concise_open, order(res_prop_concise_open[[varnames_short[[i]][1]]])),]
  test = rollapply(res_prop_ordered_open, width = 2000, FUN = mean)
  test = data.frame(test)
  for (j in 1:length(varnames_short))
  {
    if (j != i)
    {
      filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_open_cells.pdf')
      pdf( filename, width=4.5, height=4, useDingbats=FALSE )
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      plot(test[[varnames_short[[i]][1]]],
           test[[varnames_short[[j]][1]]],
           xlab = varnames_long[i],
           ylab = varnames_long[j],
           type='l')
      graphics.off()
    }
  }
}

#test_quantile = rollapply(cbind(res_prop_all_ordered$wcnSC,res_prop_all_ordered$volume, res_prop_all_ordered$rsa), width = 1000, FUN = quantile)
#test = data.frame(test_quantile)
#View(test)


plot(test$vsphericity,test$seqent, type='l')
plot(test$vsphericity,test$ddgent, type='l')
plot(test$vsphericity,test$vvolume, type='l')
plot(test$vsphericity,test$vedge, type='l')
plot(test$vsphericity,test$wcnSC, type='l')
plot(test$vsphericity,test$rsa, type='l')
plot(test$vsphericity,test$bfSC, type='l')
plot(test$vsphericity,test$hbe, type='l')
plot(test$vsphericity,test$hpshh, type='l')
plot(test$vsphericity,test$veccentricity, type='l')
  
length(res_prop_all_ordered$wcnSC)
length(vol_mean$vol_mean)
