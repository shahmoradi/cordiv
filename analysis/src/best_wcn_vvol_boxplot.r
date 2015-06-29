# This R script combines parts of the two scripts best_wcn_select_variables.r & best_vvolume_select_variables.r for the purpose of generating one figure combined for the manuscript.

# This code requires resulting data from the two aforementioned scripts.
# Last updated by Amir Shahmoradi, Wednesday 2:38 PM, June 24 2015, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

crd_list = c('SC','AA','CB','CA','C','N','O') # list of representative coordinates used to generate Voronoi cells.

# Create box plots all in one figure

  wcn_scors_all_pdbs = read.csv( "../tables/best_wcn/select_variables/wcn_scors_all_pdbs.csv", header = TRUE )
  wcn_scors_all_pdbs$variable = factor(wcn_scors_all_pdbs$variable)
  wcn_scors_all_pdbs$wcn = factor(wcn_scors_all_pdbs$wcn)

  vvolume_scors_all_pdbs = read.csv( "../tables/best_vvolume/select_variables/vvolume_scors_all_pdbs.csv", header = TRUE )
  vvolume_scors_all_pdbs$variable = factor(vvolume_scors_all_pdbs$variable)
  vvolume_scors_all_pdbs$vvolume = factor(vvolume_scors_all_pdbs$vvolume)
  
  counter = 0
  filename = paste0('../figures/best_wcn_vvol_boxplot.pdf')
  pdf( filename, width=15, height=8, useDingbats=FALSE )
  split.screen(c(2,2))
  
  screen(1)
    temp_data_wcn_long = wcn_scors_all_pdbs[wcn_scors_all_pdbs$variable == 'r4s_JC',c('pdb','wcn','value')]
    temp_data_wcn = reshape(temp_data_wcn_long, timevar = 'wcn', idvar = 'pdb', direction = 'wide')
    temp_data_wcn = subset (temp_data_wcn, select = -c(pdb))
    colnames(temp_data_wcn) = substr(levels(wcn_scors_all_pdbs$wcn),start=4,stop=5)
    temp_data_wcn = temp_data_wcn[crd_list]
    par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    boxplot(  temp_data_wcn
            , xlab = 'Representative Weighted Contact Number (WCN)'
            , ylab = paste0('Correlation with ER')
            , cex.axis = 1.3
            , cex.lab = 1.3
            , ylim = c(-0.1,-1.0)
            , outline=FALSE    # Remove outliers from the plot
            )
    mtext('A', side = 3, at=-0.35, font=2, cex=1.2)

  screen(2)
    temp_data_vvolume_long = vvolume_scors_all_pdbs[vvolume_scors_all_pdbs$variable == 'r4s_JC',c('pdb','vvolume','value')]
    temp_data_vvolume = reshape(temp_data_vvolume_long, timevar = 'vvolume', idvar = 'pdb', direction = 'wide')
    temp_data_vvolume = subset (temp_data_vvolume, select = -c(pdb))
    colnames(temp_data_vvolume) = substr(levels(vvolume_scors_all_pdbs$vvolume),start=8,stop=9)
    temp_data_vvolume = temp_data_vvolume[crd_list]
    par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    boxplot(  temp_data_vvolume
            , xlab = 'Representative Voronoi Cell Volume'
            , ylab = paste0('Correlation with ER')
            , cex.axis = 1.3
            , cex.lab = 1.3
            , ylim = c(0.1,1.)
            , outline=FALSE    # Remove outliers from the plot
            )
    mtext('B', side = 3, at=-0.35, font=2, cex=1.2)
  
  screen(3)
    temp_data_wcn_long = wcn_scors_all_pdbs[wcn_scors_all_pdbs$variable == 'rsa',c('pdb','wcn','value')]
    temp_data_wcn = reshape(temp_data_wcn_long, timevar = 'wcn', idvar = 'pdb', direction = 'wide')
    temp_data_wcn = subset (temp_data_wcn, select = -c(pdb))
    colnames(temp_data_wcn) = substr(levels(wcn_scors_all_pdbs$wcn),start=4,stop=5)
    temp_data_wcn = temp_data_wcn[crd_list]
    par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    boxplot(  temp_data_wcn
            , xlab = 'Representative Weighted Contact Number (WCN)'
            , ylab = paste0('Correlation with RSA')
            , cex.axis = 1.3
            , cex.lab = 1.3
            , ylim = c(-0.1,-1.0)
            , outline=FALSE    # Remove outliers from the plot
            )
    mtext('C', side = 3, at=-0.35, font=2, cex=1.2)
  
  screen(4)
    temp_data_vvolume_long = vvolume_scors_all_pdbs[vvolume_scors_all_pdbs$variable == 'rsa',c('pdb','vvolume','value')]
    temp_data_vvolume = reshape(temp_data_vvolume_long, timevar = 'vvolume', idvar = 'pdb', direction = 'wide')
    temp_data_vvolume = subset (temp_data_vvolume, select = -c(pdb))
    colnames(temp_data_vvolume) = substr(levels(vvolume_scors_all_pdbs$vvolume),start=8,stop=9)
    temp_data_vvolume = temp_data_vvolume[crd_list]
    par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    boxplot(  temp_data_vvolume
            , xlab = 'Representative Voronoi Cell Volume'
            , ylab = paste0('Correlation with RSA')
            , cex.axis = 1.3
            , cex.lab = 1.3
            , ylim = c(0.1,1.)
            , outline=FALSE    # Remove outliers from the plot
            )
    mtext('D', side = 3, at=-0.35, font=2, cex=1.2)

graphics.off()
