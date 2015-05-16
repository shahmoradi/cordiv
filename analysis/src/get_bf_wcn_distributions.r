# Amir Shahmoradi, Friday 8:26 PM, May 11 2015, Wilke Lab, ICMB, UT Austin

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
res_prop_wcnSC_bfSC = data.frame( pdb   = res_prop_wcn_bf$pdb
                                , wcnSC = res_prop_wcn_bf$wcnSC
                                , bfSC  = res_prop_wcn_bf$bfSC
                                )
res_prop_wcnSC_bfSC$pdb = factor(res_prop_wcnSC_bfSC$pdb)

# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Side-Chain Weighted Contact Number' , 'Average Side-Chain B factor')


varnames_short = colnames(res_prop_wcnSC_bfSC)

counter = 0
for (pdb in levels(res_prop_wcnSC_bfSC$pdb))
{
  counter = counter + 1
  cat('generating Figures for pdb ', counter, pdb, '\n')
  temp = res_prop_wcnSC_bfSC[res_prop_wcnSC_bfSC$pdb==pdb,]
  
  # B factor plots
  hist.bfSC  = density(temp$bfSC)
  filename = paste0('../figures/wcn_bf_hist/bf/hist_bfSC_',pdb,'.png')
  png( filename, width=400, height=400 )
  par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot( hist.bfSC$x
      , hist.bfSC$y
      , col = 'red'
      , type = 'l'
      , lwd  = 3 
      , xlab = expression( 'B factor [ Å'^2*' ]')
      , ylab = 'Normalized Frequency'
      )
  legend( 'topright' , paste0('PDB ',pdb) , bty = 'n' )
  graphics.off()
  
  # WCN plots
  hist.wcnSC = density(temp$wcnSC)
  filename = paste0('../figures/wcn_bf_hist/wcn/hist_bfSC_',pdb,'.png')
  png( filename, width=400, height=400 )
  par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot( hist.wcnSC$x
        , hist.wcnSC$y
        , col = 'red'
        , type = 'l'
        , lwd  = 3 
        , xlab = expression( 'WCN [ Å'^-2*' ]')
        , ylab = 'Normalized Frequency'
  )
  legend( 'topright' , paste0('PDB ',pdb) , bty = 'n' )
  graphics.off()
}

