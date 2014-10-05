# This code is aimed at finding which definition of Voronoi cell volumes (vorvol) can best represent a residue. This is done by comparing the correlations of different vorvols with other residue variables.
# "selected_variables" in the name of this code refers to the fact that I am comparing the performance of different vorvols to only a select number of important residue properties, such as ASA, RSA, seqent, ddGent, and I am ignoring the rest which are generally definition dependent, such as different definitions of B factors and WCN based on the set of atom coordinates used.

# Last updated by Amir Shahmoradi, Friday 1:23 PM, October 4 2014, Wilke Lab, ICMB, UT Austin

# Input data:
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out
#               ../../properties/res_prop_voronoiAA.out
#               ../../properties/res_prop_voronoiCA.out
#               ../../properties/res_prop_voronoiSC.out

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

vorvol_scors_all_pdbs = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
vorvol_list = c('vSC','vAA','vCA','mvSC','mvAA','mvCA','fvSC','fvAA','fvCA')
vorvol_list_long = c('VSCvolume','VAAvolume','VCAvolume','VSCmodified_volume','VAAmodified_volume','VCAmodified_volume','VSCfree_volume','VAAfree_volume','VCAfree_volume')
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,]  # c('hpskd','hpsww','hpshh')]
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  #pdb_bf     = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb]
  #pdb_wcn    = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb]
  #pdb_voroAA = res_prop_voroAA[res_prop_voroAA$pdb==pdb, ]
  #pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  #pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  pdb_voro = cbind(res_prop_voroAA[res_prop_voroAA$pdb==pdb, c('VAAvolume','VAAmodified_volume','VAAfree_volume')],
                   res_prop_voroCA[res_prop_voroCA$pdb==pdb, c('VCAvolume','VCAmodified_volume','VCAfree_volume')],
                   res_prop_voroSC[res_prop_voroSC$pdb==pdb, c('VSCvolume','VSCmodified_volume','VSCfree_volume')])
  pdb_voro = pdb_voro[, vorvol_list_long]
  colnames(pdb_voro) = vorvol_list
  
  pdb_temp = cbind( subset(pdb_elj, select  = c(seqent,ddgent)),
                    subset(pdb_hps, select  = c(hpshh)),
                    subset(pdb_dssp, select = c(asa,rsa,hbe_mean))
                    #subset(pdb_wcn, select = -c(pdb,resnam,resnum,bfSC,bfAA,bfN,bfCA,bfC,bfO,bfCB)),
                    #subset(pdb_voroAA, select = c(VAAvolume,VAAmodified_volume,VAAfree_volume)),
                    #subset(pdb_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change)),
                    #subset(pdb_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change))
                    )
  
  pdb_vorvol_long = reshape(pdb_voro, ids = rownames(pdb_voro), varying = colnames(pdb_voro), v.names = 'value', timevar = 'variable', times = colnames(pdb_voro), direction = 'long')
  pdb_vorvol_long$variable = factor(pdb_vorvol_long$variable)
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  for (vorvol in levels(pdb_vorvol_long$variable))
  {
    vorvol_data = pdb_vorvol_long[pdb_vorvol_long$variable == vorvol,]
    for (variable in levels(pdb_long$variable))
    {
      variable_data = pdb_long[pdb_long$variable == variable,]
      x = cor.test( vorvol_data$value, variable_data$value, method='spearman', na.action="na.omit" )
      r = x$estimate
      p = x$p.value
      
      row = data.frame(pdb = pdb, vorvol = vorvol, variable = variable, value = r)
      vorvol_scors_all_pdbs = rbind(vorvol_scors_all_pdbs,row)
    }
  }
}
write.csv( vorvol_scors_all_pdbs, "../tables/best_vorvol/selected_variables/vorvol_scors_all_pdbs.csv", row.names=F )



# Now summarize all Spearman correlations over the entire dataset

vorvol_scors_all_pdbs = read.csv( "../tables/best_vorvol/selected_variables/vorvol_scors_all_pdbs.csv", header = TRUE )
vorvol_scors_all_pdbs$variable = factor(vorvol_scors_all_pdbs$variable)
vorvol_scors_all_pdbs$vorvol = factor(vorvol_scors_all_pdbs$vorvol)


for (vorvol in levels(vorvol_scors_all_pdbs$vorvol))
{
  temp_data_vorvol = vorvol_scors_all_pdbs[vorvol_scors_all_pdbs$vorvol == vorvol,]
  temp_data_vorvol$variable = factor(temp_data_vorvol$variable)
  
  vorvol_scors_summary = data.frame()
  for (variable in levels(temp_data_vorvol$variable))
  {
    temp_data = temp_data_vorvol[temp_data_vorvol$variable == variable,]
    if (length(temp_data$pdb) != 213)
    {
      stop ( 'something is fishy here!' )
    }
    
    x   = quantile(temp_data$value, probs = c(0,0.25,0.5,0.75, 1.))
    row = data.frame(variable = variable,
                     mean     = mean(temp_data$value),
                     median   = x[3],
                     sd       = sd(temp_data$value),
                     min      = x[1],
                     q1       = x[2],
                     q3       = x[4],
                     max      = x[5]
                    )
    vorvol_scors_summary = rbind(vorvol_scors_summary,row)
  }
  row.names(vorvol_scors_summary) = c()
  filename = paste0("../tables/best_vorvol/selected_variables/",vorvol,'_scors_summary.csv')
  write.csv(vorvol_scors_summary, filename, row.names=F )
  cat (vorvol, filename, '\n')
}

# Now compare the three vorvol scors with each other:

for (i in 1:(length(vorvol_list)-1))
{
  filename = paste0("../tables/best_vorvol/selected_variables/",vorvol_list[[i]][1],'_scors_summary.csv')
  vorvol_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(vorvol_list))
  {
    filename = paste0("../tables/best_vorvol/selected_variables/",vorvol_list[[j]][1],'_scors_summary.csv')
    vorvol_scors_summary2 = read.csv( filename, header = T )
    difference = data.frame(variable      = vorvol_scors_summary1$variable,
                            mean_diff     = abs(vorvol_scors_summary1$mean)   - abs(vorvol_scors_summary2$mean),
                            median_diff   = abs(vorvol_scors_summary1$median) - abs(vorvol_scors_summary2$median),
                            sd_diff       = abs(vorvol_scors_summary1$sd)     - abs(vorvol_scors_summary2$sd),
                            min_diff      = abs(vorvol_scors_summary1$min)    - abs(vorvol_scors_summary2$min),
                            q1_diff       = abs(vorvol_scors_summary1$q1)     - abs(vorvol_scors_summary2$q1),
                            q3_diff       = abs(vorvol_scors_summary1$q3)     - abs(vorvol_scors_summary2$q3),
                            max_diff      = abs(vorvol_scors_summary1$max)    - abs(vorvol_scors_summary2$max)
                            )
    row.names(difference) = c()
    filename = paste0('../tables/best_vorvol/selected_variables/diff_',vorvol_list[[i]][1],'_',vorvol_list[[j]][1],'.csv')
    cat (filename, '\n')
    write.csv( difference, filename, row.names=F )
  }
}

# Make plots of ABS(vorvol) vs. ABS(vorvol) for comparison:

x = -2:2
colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
labels = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')

for (i in 1:(length(vorvol_list)-1))
{
  filename = paste0("../tables/best_vorvol/selected_variables/",vorvol_list[[i]][1],'_scors_summary.csv')
  vorvol_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(vorvol_list))
  {
    filename = paste0("../tables/best_vorvol/selected_variables/",vorvol_list[[j]][1],'_scors_summary.csv')
    vorvol_scors_summary2 = read.csv( filename, header = T )
    filename = paste0('../figures/best_vorvol/selected_variables/abs(',vorvol_list[[i]][1],'_',vorvol_list[[j]][1],').pdf')
    cat (filename, '\n')
    pdf( filename, width=4.5, height=4, useDingbats=FALSE )
    par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    plot(-999,
         #xaxt='n',yaxt='n',bty='n',pch='',
    #plot(abs(vorvol_scors_summary1$median),
    #     abs(vorvol_scors_summary2$median),
         xlab = paste0( 'absolute median correlation with ',vorvol_list[[i]][1] ),
         ylab = paste0( 'absolute median correlation with ',vorvol_list[[j]][1] ),
         xlim = c(0,1.0),
         ylim = c(0,1.0)
         )
    points( abs(vorvol_scors_summary1$median),
            abs(vorvol_scors_summary2$median), pch=19, col=colors[1:length(vorvol_scors_summary2$median)])
    lines( x, x, col='red' )
    legend( -0.02, 1.05, labels[1:7], pch=19, col=colors[1:7], bty='n', cex=0.9)
    #legend( 0.7, 0.1, labels[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
    graphics.off()
  }
}

# Now create box plots

vorvol_scors_all_pdbs = read.csv( "../tables/best_vorvol/selected_variables/vorvol_scors_all_pdbs.csv", header = TRUE )
vorvol_scors_all_pdbs$variable = factor(vorvol_scors_all_pdbs$variable)
vorvol_scors_all_pdbs$vorvol = factor(vorvol_scors_all_pdbs$vorvol)

variable_list = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')
counter = 0
for (variable in levels(vorvol_scors_all_pdbs$variable))
{ 
  counter = counter + 1
  temp_data_variable_long = vorvol_scors_all_pdbs[vorvol_scors_all_pdbs$variable == variable,c('pdb','vorvol','value')]
  temp_data_variable = reshape(temp_data_variable_long, timevar = 'vorvol', idvar = 'pdb', direction = 'wide')
  temp_data_variable = subset (temp_data_variable, select = -c(pdb))
  colnames(temp_data_variable) = levels(vorvol_scors_all_pdbs$vorvol)
  temp_data_variable = temp_data_variable[vorvol_list]
  filename = paste0('../figures/best_vorvol/selected_variables/boxplot_',variable,'.pdf')
  pdf( filename, width=6, height=4, useDingbats=FALSE )
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  boxplot(temp_data_variable,
          xlab = 'representative Voronoi cell volume (vorvol)',
          ylab = paste0('Spearman correlation with ',variable_list[[counter]][1])
  )
  graphics.off()
}

