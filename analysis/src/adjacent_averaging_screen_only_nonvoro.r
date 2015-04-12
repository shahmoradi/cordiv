# Amir Shahmoradi, Saturday 3:07 PM, April 11 2015, Wilke Lab, ICMB, UT Austin
# This code is very similar to the original code adjacent_averaging.r with the only difference being that here the figures are generated as multiple plots for the a single residue variable versus all others, in a single pdf figure.
# This code only makes screen plots of 3*2 size

#install.packages('zoo')
#install.packages('ggplot2')
#install.packages("fields")
library('ggplot2')
library('fields') # used for function fudgeit
library('zoo')


# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# source('get_res_data.r')

# The following is short version of the most useful residue properties that I have found so far.
# Now since the three Voronoi quantities vnvortices, vnedges and vnfaces happen to be exactly the same as seen in the cormat calculated above, I am going to remove them from the data set, in addition to resvol which is not as informative.
res_prop_concise = data.frame( zr4s_JC       = res_prop_jec$zr4s_JC
                             , seqent        = res_prop_elj$seqent
                             , ddgent        = res_prop_jec_ddg$rate.ddg.foldx
                             , rsa           = res_prop_dssp$rsa
                             , wcnSC         = res_prop_wcn_bf$wcnSC
                             , bfSC          = res_prop_wcn_bf$bfSC
                             , hbe           = res_prop_dssp$hbe
                             )

res_prop_concise_closed = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff == 0,]
res_prop_concise_open = res_prop_concise[res_prop_voroSC$VSCvolume_change_diff != 0,]

#qplot(wcnSC, zr4s_JC, data = res_prop_concise, geom = c("smooth"))
#qplot(rsa, zr4s_JC, data = res_prop_concise, geom = c("smooth"))
#qplot(varea, zr4s_JC, data = res_prop_concise, geom = c("smooth"))
#qplot(ddgent, zr4s_JC, data = res_prop_concise, geom = c("smooth"))
#qplot(log10(bfSC), zr4s_JC, data = res_prop_concise, geom = c("point"))


# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Evolutionary Rates (r4sJC)' , 'Sequence Entropy' , 'ddG Rate (Stability upon Substitution)' , 'Relative Solvent Accessibility (RSA)' ,
				  'Hydrophobicity Scale' , 'Side-Chain Weighted Contact Number' , 'Average Side-Chain B factor' , 'Hydrogen Bond Energy (HBE)')


varnames_short = colnames(res_prop_concise)


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


fudgeit <- function(){
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  image.plot(xm,ym,z, col = colramp(256)
             , legend.only = T
             , add = F
  )
}


for (i in 1:length(varnames_short))
{
  #i = 1
  cat('generating graph # ', i, '\n')
  filename = paste0('../figures/adjacent_averaging_screen/only_nonvoro/',varnames_short[i],'_nonvoro.png')
  #pdf( filename, width=13.5, height=16, useDingbats=FALSE )  
  png( filename, width=740, height=800 )  
  split.screen(c(3,2))
  
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
      par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 ) #, mar = c(5,4,4,5) + .1 )
      
      smoothScatter(temp[[varnames_short[[j]][1]]],
                    temp[[varnames_short[[i]][1]]],
                    xlab = varnames_long[j],
                    ylab = varnames_long[i],
                    nrpoints=0,
                    nbin=500
                    #,postPlotHook = fudgeit
                    ) 
      lines(temp[[varnames_short[[j]][1]]],
           temp[[varnames_short[[i]][1]]],
           col = 'red',
           #type='l'
           lwd=3
           )
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
