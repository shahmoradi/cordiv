# Amir Shahmoradi, Sunday 7:42 PM, Oct 5 2014, Wilke Lab, ICMB, UT Austin
# This code is very similar to the original code adjacent_averaging.r with the only difference being that here the figures are generated as multiple plots for the a single residue variable versus all others, in a single pdf figure.

# install.packages('zoo')
library('zoo')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# The following is short version of the most useful residue properties that I have found so far.
# Now since the three Voronoi quantities vnvortices, vnedges and vnfaces happen to be exactly the same as seen in the cormat calculated above, I am going to remove them from the data set, in addition to resvol which is not as informative.
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
                              vedge         = log10(res_prop_voroSC$VSCedge_length_total),
                              varea         = log10(res_prop_voroSC$VSCarea),
                              vvolume       = log10(res_prop_voroSC$VSCvolume),
                              veccentricity = res_prop_voroSC$VSCeccentricity,
                              vsphericity   = res_prop_voroSC$VSCsphericity
                              #mvsphericity  = res_prop_voroSC$VSCsphericity
                              )

res_prop_concise_closed = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff == 0,]
res_prop_concise_open = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff != 0,]

# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Sequence Entropy (seqent)' , 'ddG Entropy (ddGent)' , 'Relative Solvent Accessibility (RSA)' , 'Hydrogen Bond Energy (HBE)' ,
                  'Hydrophobicity Scale (HPS)' , 'Side-Chain Contact Number (wcnSC)' ,
                  'Average Side-Chain Bfactor (bfSC)' , 'Voronoi Cell Faces' , 'log10 ( Voronoi Cell Edge length )' , 'log10 ( Voronoi Cell Surface Area )' ,
                  'log10 ( Voronoi Cell Volume )' , 'Voronoi Cell Eccentricity' , 'Voronoi Cell Sphericity')

varnames_short = colnames(res_prop_concise)

voronoi_colnames = c( 'vnfaces' , 'vedge' , 'varea' , 'vvolume' , 'veccentricity' , 'vsphericity')

# In order to reduce runtime, I am going to save the ordered data frames in a list once, which will be then used in the loops wherever needed.
list_temp = list()
list_temp_closed = list()
list_temp_open   = list()
for (i in 1:length(varnames_short))
{
  cat('generating list element ', i, '\n')
  res_prop_ordered = res_prop_concise[with(res_prop_concise, order(res_prop_concise[[varnames_short[[i]][1]]])),]
  temp = rollapply(res_prop_ordered, width = 3000, FUN = mean)
  temp = data.frame(temp)
  list_temp[[i]] = temp
  
  res_prop_ordered_closed = res_prop_concise_closed[with(res_prop_concise_closed, order(res_prop_concise_closed[[varnames_short[[i]][1]]])),]
  temp_closed = rollapply(res_prop_ordered_closed, width = 3000, FUN = mean)
  temp_closed = data.frame(temp_closed)
  list_temp_closed[[i]] = temp_closed
  
  res_prop_ordered_open = res_prop_concise_open[with(res_prop_concise_open, order(res_prop_concise_open[[varnames_short[[i]][1]]])),]
  temp_open = rollapply(res_prop_ordered_open, width = 3000, FUN = mean)
  temp_open = data.frame(temp_open)
  list_temp_open[[i]] = temp_open
}


for (i in 1:length(varnames_short))
{
  #i = 1
  cat('generating graph # ', i, '\n')
  filename = paste0('../figures/adjacent_averaging_screen/voro_log/',varnames_short[i],'.pdf')
  pdf( filename, width=13.5, height=16, useDingbats=FALSE )  
  split.screen(c(4,3))
  
  counter = 0
  for (j in 1:length(varnames_short))
  {
    
    if (j != i)
    {
      # Now call the elements of the lists generated above.
      temp = as.data.frame(list_temp[[j]])
      temp_closed = as.data.frame(list_temp_closed[[j]])
      temp_open = as.data.frame(list_temp_open[[j]])
      
      counter = counter + 1
      screen(counter)
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      if ((varnames_short[[i]][1] %in% voronoi_colnames) | (varnames_short[[j]][1] %in% voronoi_colnames))
      {
        # Plot all closed, open and all cells in one.
        x_range = range( min(temp[[varnames_short[[j]][1]]], temp_closed[[varnames_short[[j]][1]]], temp_open[[varnames_short[[j]][1]]]),
                         max(temp[[varnames_short[[j]][1]]], temp_closed[[varnames_short[[j]][1]]], temp_open[[varnames_short[[j]][1]]]) )
        y_range = range( min(temp[[varnames_short[[i]][1]]], temp_closed[[varnames_short[[i]][1]]], temp_open[[varnames_short[[i]][1]]]),
                         max(temp[[varnames_short[[i]][1]]], temp_closed[[varnames_short[[i]][1]]], temp_open[[varnames_short[[i]][1]]]) )
        plot(temp[[varnames_short[[j]][1]]],
             temp[[varnames_short[[i]][1]]],
             xlab = varnames_long[j],
             ylab = varnames_long[i],
             xlim = x_range,
             ylim = y_range,
             col = 'black',
             type='l'
             )
        lines(temp_closed[[varnames_short[[j]][1]]],
              temp_closed[[varnames_short[[i]][1]]],
              col = 'red')
        lines(temp_open[[varnames_short[[j]][1]]],
              temp_open[[varnames_short[[i]][1]]],
              col = 'blue')
      }
      else
      {
        plot(temp[[varnames_short[[j]][1]]],
             temp[[varnames_short[[i]][1]]],
             xlab = varnames_long[j],
             ylab = varnames_long[i],
             type='l')
      }
    }
  }
  close.screen(all = TRUE)
  graphics.off()
}

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
