# Amir Shahmoradi, Wednesday 3:23 PM, Sep 24 2014, Wilke Lab, ICMB, UT Austin

# install.packages('zoo')
library('zoo')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

source('get_res_data.r')

# The following is short version of the most useful residue properties that I have found so far.
# Now since the three Voronoi quantities vnvortices, vnedges and vnfaces happen to be exactly the same as seen in the cormat calculated above, I am going to remove them from the data set, in addition to resvol which is not as informative.
res_prop_concise = data.frame(zr4s_JC       = res_prop_jec$zr4s_JC,
                              seqent        = res_prop_elj$seqent,
                              ddgent        = res_prop_elj$ddgent,
                              rsa           = res_prop_dssp$rsa,
                              hbe           = res_prop_dssp$hbe,
                              hpshh         = res_prop_hps$hpshh,
                              wcnSC         = res_prop_wcn_bf$wcnSC,
                              bfSC          = res_prop_wcn_bf$bfSC,
                              #resvol        = res_prop_voroSC$resvol,
                              #vnvertices    = res_prop_voroSC$VSCnvertices,
                              #vnedges       = res_prop_voroSC$VSCnedges,
                              vnfaces       = res_prop_voroSC$VSCnfaces,
                              #vedge         = log10(res_prop_voroSC$VSCedge_length_total),
                              varea         = log10(res_prop_voroSC$VSCarea),
                              vvolume       = log10(res_prop_voroSC$VSCvolume),
                              veccentricity = res_prop_voroSC$VSCeccentricity,
                              vsphericity   = res_prop_voroSC$VSCsphericity,
                              mvsphericity  = res_prop_voroSC$VSCsphericity
                              )
# Now use the negative value of sphericity for those voronoi cells that are open.
res_prop_concise$mvsphericity[res_prop_voroSC$VSCvolume_change_diff != 0] = -res_prop_concise$vsphericity[res_prop_voroSC$VSCvolume_change_diff != 0]

cormat = cor(res_prop_concise, method='spearman')
write.csv( cormat, "../tables/res_prop_cormat.csv", row.names=T )

res_prop_concise_closed = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff == 0,]
cormat_closed = cor(res_prop_concise_closed, method='spearman')
write.csv( cormat_closed, "../tables/res_prop_cormat_closed_Vcells.csv", row.names=T )
# Now write out the difference of the two cormats for all cells and only closed cells.
cormat_diff = abs(cormat) - abs(cormat_closed)
write.csv( cormat_diff, "../tables/res_prop_cormat_all_closed_diff.csv", row.names=T )

res_prop_concise_open = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff != 0,]
cormat_open = cor(res_prop_concise_open, method='spearman')
write.csv( cormat_open, "../tables/res_prop_cormat_open_Vcells.csv", row.names=T )
# Now write out the difference of the two cormats for all cells and only closed cells.
cormat_diff = abs(cormat_open) - abs(cormat_closed)
write.csv( cormat_diff, "../tables/res_prop_cormat_open_closed_diff.csv", row.names=T )

# Now remove the modified sphericity variable from data, in order to generate plots.
res_prop_concise = subset( res_prop_concise, select = -c(mvsphericity) )
res_prop_concise_closed = subset( res_prop_concise_closed, select = -c(mvsphericity) )
res_prop_concise_open = subset( res_prop_concise_open, select = -c(mvsphericity) )


# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Evolutionary Rates (r4sJC)' , 'Sequence Entropy (seqent)' , 'ddG Entropy (ddGent)' , 'Relative Solvent Accessibility (RSA)' , 'Hydrogen Bond Energy (HBE)' ,
                  'Hydrophobicity Scale (HPS)' , 'Side-Chain Contact Number (wcnSC)' ,
                  'Average Side-Chain Bfactor (bfSC)' , 'Voronoi Cell Faces' , 'log10 ( Voronoi Cell Surface Area )' ,
                  'log10 ( Voronoi Cell Volume )' , 'Voronoi Cell Eccentricity' , 'Voronoi Cell Sphericity')

varnames_short = colnames(res_prop_concise)

voronoi_colnames = c( 'vnfaces' , 'varea' , 'vvolume' , 'veccentricity' , 'vsphericity')

for (i in 1:length(varnames_short))
{
  cat(i, varnames_short[i], '\n')
  res_prop_ordered = res_prop_concise[with(res_prop_concise, order(res_prop_concise[[varnames_short[[i]][1]]])),]
  temp = rollapply(res_prop_ordered, width = 3000, FUN = mean)
  temp = data.frame(temp)
  
  res_prop_ordered_closed = res_prop_concise_closed[with(res_prop_concise_closed, order(res_prop_concise_closed[[varnames_short[[i]][1]]])),]
  temp_closed = rollapply(res_prop_ordered_closed, width = 3000, FUN = mean)
  temp_closed = data.frame(temp_closed)
  
  res_prop_ordered_open = res_prop_concise_open[with(res_prop_concise_open, order(res_prop_concise_open[[varnames_short[[i]][1]]])),]
  temp_open = rollapply(res_prop_ordered_open, width = 3000, FUN = mean)
  temp_open = data.frame(temp_open)
  
  for (j in 1:length(varnames_short))
  {
    if (j != i)
    {

      filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_all_cells.pdf')
      pdf( filename, width=5.625, height=5, useDingbats=FALSE )
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      plot(temp[[varnames_short[[i]][1]]],
           temp[[varnames_short[[j]][1]]],
           xlab = varnames_long[i],
           ylab = varnames_long[j],
           type='l')
      graphics.off()
      
      if ((varnames_short[[i]][1] %in% voronoi_colnames) | (varnames_short[[j]][1] %in% voronoi_colnames))
      {
        filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_closed_cells.pdf')
        pdf( filename, width=5.625, height=5, useDingbats=FALSE )
        par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
        plot(temp_closed[[varnames_short[[i]][1]]],
             temp_closed[[varnames_short[[j]][1]]],
             xlab = varnames_long[i],
             ylab = varnames_long[j],
             type='l')
        graphics.off()
        
        filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_open_cells.pdf')
        pdf( filename, width=5.625, height=5, useDingbats=FALSE )
        par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
        plot(temp_open[[varnames_short[[i]][1]]],
             temp_open[[varnames_short[[j]][1]]],
             xlab = varnames_long[i],
             ylab = varnames_long[j],
             type='l')
        graphics.off()
        
        # Now plot all three in one.
        filename = paste0('../figures/adjacent_averaging/',varnames_short[i],'_',varnames_short[j],'_all_in_one.pdf')
        x_range = range( min(temp[[varnames_short[[i]][1]]], temp_closed[[varnames_short[[i]][1]]], temp_open[[varnames_short[[i]][1]]]),
                         max(temp[[varnames_short[[i]][1]]], temp_closed[[varnames_short[[i]][1]]], temp_open[[varnames_short[[i]][1]]]) )
        y_range = range( min(temp[[varnames_short[[j]][1]]], temp_closed[[varnames_short[[j]][1]]], temp_open[[varnames_short[[j]][1]]]),
                         max(temp[[varnames_short[[j]][1]]], temp_closed[[varnames_short[[j]][1]]], temp_open[[varnames_short[[j]][1]]]) )
        pdf( filename, width=5.625, height=5, useDingbats=FALSE )
        par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
        plot(temp[[varnames_short[[i]][1]]],
             temp[[varnames_short[[j]][1]]],
             xlab = varnames_long[i],
             ylab = varnames_long[j],
             xlim = x_range,
             ylim = y_range,
             col = 'black',
             type='l'
             )
        lines(temp_closed[[varnames_short[[i]][1]]],
              temp_closed[[varnames_short[[j]][1]]],
              col = 'red')
        lines(temp_open[[varnames_short[[i]][1]]],
              temp_open[[varnames_short[[j]][1]]],
              col = 'blue')
        graphics.off()
      }
    }
  }
}

# Update: Feb 25 2015, Amir
# Now let's make a plot of r4sJC vs. inverse of wcnSC, that is 1/wcnSC, to see if it is linear or not.

res_prop_ordered = res_prop_concise[with(res_prop_concise, order(res_prop_concise$wcnSC)),]
temp = rollapply(res_prop_ordered, width = 3000, FUN = mean)
temp = data.frame(temp)
filename = paste0('../figures/adjacent_averaging/wcnSCinv_zr4sJC_all_cells.pdf')
pdf( filename, width=5.625, height=5, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(1./temp$wcnSC,
     temp$zr4s_JC,
     xlab = "1 / wcnSC",
     ylab = "Evolutionary Rates (zr4sJC)",
     type='l')
graphics.off()

# Update: Feb 26 2015, Amir
# Now let's make a plot of r4sJC vs. inverse of Voronoi Volume (volinv), that is 1/volume, to see if it is linear or not.

res_prop_ordered = res_prop_concise[with(res_prop_concise, order(res_prop_concise$vvolume)),]
temp = rollapply(res_prop_ordered, width = 3000, FUN = mean)
temp = data.frame(temp)
filename = paste0('../figures/adjacent_averaging/vvolumeSCinv_zr4sJC_all_cells.pdf')
pdf( filename, width=5.625, height=5, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(1./temp$vvolume,
     temp$zr4s_JC,
     xlab = "1 / Voronoi Volume",
     ylab = "Evolutionary Rates (zr4sJC)",
     type='l')
graphics.off()

res_prop_concise_closed = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff == 0,]
res_prop_closed_ordered = res_prop_concise_closed[with(res_prop_concise_closed, order(res_prop_concise_closed$vvolume)),]
temp_closed = data.frame(rollapply(res_prop_closed_ordered, width = 3000, FUN = mean))
filename = paste0('../figures/adjacent_averaging/vvolumeSCinv_zr4sJC_all_cells_closed_cells.pdf')
pdf( filename, width=5.625, height=5, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(1./temp_closed$vvolume,
     temp_closed$zr4s_JC,
     xlab = "1 / Voronoi Volume",
     ylab = "Evolutionary Rates (zr4sJC)",
     type='l')
graphics.off()

res_prop_concise_open = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff != 0,]
res_prop_open_ordered = res_prop_concise_open[with(res_prop_concise_open, order(res_prop_concise_open$vvolume)),]
temp_open = data.frame(rollapply(res_prop_open_ordered, width = 3000, FUN = mean))
filename = paste0('../figures/adjacent_averaging/vvolumeSCinv_zr4sJC_all_cells_open_cells.pdf')
pdf( filename, width=5.625, height=5, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(1./temp_open$vvolume,
     temp_open$zr4s_JC,
     xlab = "1 / Voronoi Volume",
     ylab = "Evolutionary Rates (zr4sJC)",
     type='l')
graphics.off()



#temp_quantile = rollapply(cbind(res_prop_all_ordered$wcnSC,res_prop_all_ordered$volume, res_prop_all_ordered$rsa), width = 1000, FUN = quantile)
#temp = data.frame(temp_quantile)
#View(temp)


# plot(temp$vsphericity,temp$seqent, type='l')
# plot(temp$vsphericity,temp$ddgent, type='l')
# plot(temp$vsphericity,temp$vvolume, type='l')
# plot(temp$vsphericity,temp$vedge, type='l')
# plot(temp$vsphericity,temp$wcnSC, type='l')
# plot(temp$vsphericity,temp$rsa, type='l')
# plot(temp$vsphericity,temp$bfSC, type='l')
# plot(temp$vsphericity,temp$hbe, type='l')
# plot(temp$vsphericity,temp$hpshh, type='l')
# plot(temp$vsphericity,temp$veccentricity, type='l')
  
# length(res_prop_all_ordered$wcnSC)
# length(vol_mean$vol_mean)

#volume_smooth = filter(res_prop_all_ordered$volume, ma1000)
#plot(res_prop_all_ordered$wcnSC,volume_smooth, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )
#plot(smoothed$wcnSC_mean,smoothed$vol_mean, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )
