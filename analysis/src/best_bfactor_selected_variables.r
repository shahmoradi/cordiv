# This code is aimed at finding which definition of B factor can best represent a residue. This is done by comparing the correlations of different atomic B factors or defenitions of B factor (SC and AA) with other residue variables.
# "selected_variables" in the name of this code refers to the fact that I am comparing the performance of different B factors to only a select number of important residue properties, such as ASA, RSA, seqent, ddGent, and I am ignoring the rest which are generally definition dependent, such as different definitions of WCN and voronoi Volumes and areas based on the set of atom coordinates used.

# Last updated by Amir Shahmoradi, Wednesday 5:06 PM, October 1 2014, Wilke Lab, ICMB, UT Austin

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

bf_scors_all_pdbs = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
bf_list = c('bfSC','bfAA','bfCB','bfCA','bfN','bfC','bfO')
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,]  # c('hpskd','hpsww','hpshh')]
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_bf     = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, bf_list]
  #pdb_wcn    = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb,]
  pdb_voroAA = res_prop_voroAA[res_prop_voroAA$pdb==pdb, ]
  #pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  #pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select  = c(seqent,ddgent)),
                    subset(pdb_hps, select  = c(hpshh)),
                    subset(pdb_dssp, select = c(asa,rsa,hbe_mean)),
                    subset(pdb_voroAA, select = c(resvol))
                    #subset(pdb_wcn, select = -c(pdb,resnam,resnum,bfSC,bfAA,bfN,bfCA,bfC,bfO,bfCB)),
                    #subset(pdb_voroAA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,VAAnvertices,VAAnedges,VAAvolume_change)),
                    #subset(pdb_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change)),
                    #subset(pdb_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change))
                   )
  
  pdb_bf_long = reshape(pdb_bf, ids = rownames(pdb_bf), varying = colnames(pdb_bf), v.names = 'value', timevar = 'variable', times = colnames(pdb_bf), direction = 'long')
  pdb_bf_long$variable = factor(pdb_bf_long$variable)
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  for (bf in levels(pdb_bf_long$variable))
  {
    bf_data = pdb_bf_long[pdb_bf_long$variable == bf,]
    for (variable in levels(pdb_long$variable))
    {
      variable_data = pdb_long[pdb_long$variable == variable,]
      x = cor.test( bf_data$value, variable_data$value, method='spearman', na.action="na.omit" )
      r = x$estimate
      p = x$p.value
      
      row = data.frame(pdb = pdb, bf = bf, variable = variable, value = r)
      bf_scors_all_pdbs = rbind(bf_scors_all_pdbs,row)
    }
  }
}
write.csv( bf_scors_all_pdbs, "../tables/best_bf/selected_variables/bf_scors_all_pdbs.csv", row.names=F )



# Now summarize all Spearman correlations over the entire dataset

bf_scors_all_pdbs = read.csv( "../tables/best_bf/selected_variables/bf_scors_all_pdbs.csv", header = TRUE )
bf_scors_all_pdbs$variable = factor(bf_scors_all_pdbs$variable)
bf_scors_all_pdbs$bf = factor(bf_scors_all_pdbs$bf)


for (bf in levels(bf_scors_all_pdbs$bf))
{
  temp_data_bf = bf_scors_all_pdbs[bf_scors_all_pdbs$bf == bf,]
  temp_data_bf$variable = factor(temp_data_bf$variable)
  
  bf_scors_summary = data.frame()
  for (variable in levels(temp_data_bf$variable))
  {
    temp_data = temp_data_bf[temp_data_bf$variable == variable,]
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
    bf_scors_summary = rbind(bf_scors_summary,row)
  }
  row.names(bf_scors_summary) = c()
  filename = paste0("../tables/best_bf/selected_variables/",bf,'_scors_summary.csv')
  write.csv(bf_scors_summary, filename, row.names=F )
  cat (bf, filename, '\n')
}

# Now compare the three bf scors with each other:

for (i in 1:(length(bf_list)-1))
{
  filename = paste0("../tables/best_bf/selected_variables/",bf_list[[i]][1],'_scors_summary.csv')
  bf_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(bf_list))
  {
    filename = paste0("../tables/best_bf/selected_variables/",bf_list[[j]][1],'_scors_summary.csv')
    bf_scors_summary2 = read.csv( filename, header = T )
    difference = data.frame(variable      = bf_scors_summary1$variable,
                            mean_diff     = abs(bf_scors_summary1$mean)   - abs(bf_scors_summary2$mean),
                            median_diff   = abs(bf_scors_summary1$median) - abs(bf_scors_summary2$median),
                            sd_diff       = abs(bf_scors_summary1$sd)     - abs(bf_scors_summary2$sd),
                            min_diff      = abs(bf_scors_summary1$min)    - abs(bf_scors_summary2$min),
                            q1_diff       = abs(bf_scors_summary1$q1)     - abs(bf_scors_summary2$q1),
                            q3_diff       = abs(bf_scors_summary1$q3)     - abs(bf_scors_summary2$q3),
                            max_diff      = abs(bf_scors_summary1$max)    - abs(bf_scors_summary2$max)
                            )
    row.names(difference) = c()
    filename = paste0('../tables/best_bf/selected_variables/diff_',bf_list[[i]][1],'_',bf_list[[j]][1],'.csv')
    cat (filename, '\n')
    write.csv( difference, filename, row.names=F )
  }
}

# Make plots of ABS(bf) vs. ABS(bf) for comparison:

x = -1:1
colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
labels = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')

for (i in 1:(length(bf_list)-1))
{
  filename = paste0("../tables/best_bf/selected_variables/",bf_list[[i]][1],'_scors_summary.csv')
  bf_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(bf_list))
  {
    filename = paste0("../tables/best_bf/selected_variables/",bf_list[[j]][1],'_scors_summary.csv')
    bf_scors_summary2 = read.csv( filename, header = T )
    filename = paste0('../figures/best_bf/selected_variables/abs(',bf_list[[i]][1],'_',bf_list[[j]][1],').pdf')
    cat (filename, '\n')
    pdf( filename, width=4.5, height=4, useDingbats=FALSE )
    par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    plot(-999,
         #xaxt='n',yaxt='n',bty='n',pch='',
    #plot(abs(bf_scors_summary1$median),
    #     abs(bf_scors_summary2$median),
         xlab = paste0( 'absolute median correlation with ',bf_list[[i]][1] ),
         ylab = paste0( 'absolute median correlation with ',bf_list[[j]][1] ),
         xlim = c(0,0.8),
         ylim = c(0,0.8)
         )
    points( abs(bf_scors_summary1$median),
            abs(bf_scors_summary2$median), pch=19, col=colors[1:length(bf_scors_summary2$median)])
    lines( x, x, col='red' )
    legend( -0.02, 0.85, labels[1:7], pch=19, col=colors[1:7], bty='n', cex=0.9)
    #legend( 0.7, 0.1, labels[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
    graphics.off()
  }
}

# Now create box plots

bf_scors_all_pdbs = read.csv( "../tables/best_bf/selected_variables/bf_scors_all_pdbs.csv", header = TRUE )
bf_scors_all_pdbs$variable = factor(bf_scors_all_pdbs$variable)
bf_scors_all_pdbs$bf = factor(bf_scors_all_pdbs$bf)

variable_list = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')
counter = 0
for (variable in levels(bf_scors_all_pdbs$variable))
{ 
  counter = counter + 1
  temp_data_variable_long = bf_scors_all_pdbs[bf_scors_all_pdbs$variable == variable,c('pdb','bf','value')]
  temp_data_variable = reshape(temp_data_variable_long, timevar = 'bf', idvar = 'pdb', direction = 'wide')
  temp_data_variable = subset (temp_data_variable, select = -c(pdb))
  colnames(temp_data_variable) = levels(bf_scors_all_pdbs$bf)
  temp_data_variable = temp_data_variable[bf_list]
  filename = paste0('../figures/best_bf/selected_variables/boxplot_',variable,'.pdf')
  pdf( filename, width=6, height=4, useDingbats=FALSE )
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  boxplot(temp_data_variable,
          xlab = 'representative B factor (bf)',
          ylab = paste0('Spearman correlation with ',variable_list[[counter]][1])
          )
  graphics.off()
}

# Now create box plots all in one figure

bf_scors_all_pdbs = read.csv( "../tables/best_bf/selected_variables/bf_scors_all_pdbs.csv", header = TRUE )
bf_scors_all_pdbs$variable = factor(bf_scors_all_pdbs$variable)
bf_scors_all_pdbs$bf = factor(bf_scors_all_pdbs$bf)

variable_list = c('Sequence Entropy','ddG Entropy','RSA','ASA','Hydrophobicity','H-bond Energy')
variable_names = c('seqent','ddgent','rsa','asa','hpshh','hbe_mean')
bf_list = c('bfSC','bfAA','bfCB','bfCA','bfN','bfC','bfO')
counter = 0
filename = paste0('../figures/best_bf/selected_variables/boxplot_bf_all_in_one.pdf')
pdf( filename, width=15, height=12, useDingbats=FALSE )
split.screen(c(3,2))
for (variable in variable_names)
{ 
  counter = counter + 1
  screen(counter)
  temp_data_variable_long = bf_scors_all_pdbs[bf_scors_all_pdbs$variable == variable,c('pdb','bf','value')]
  temp_data_variable = reshape(temp_data_variable_long, timevar = 'bf', idvar = 'pdb', direction = 'wide')
  temp_data_variable = subset (temp_data_variable, select = -c(pdb))
  colnames(temp_data_variable) = levels(bf_scors_all_pdbs$bf)
  temp_data_variable = temp_data_variable[bf_list]
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  boxplot(temp_data_variable,
          xlab = 'representative B factor (bf)',
          ylab = paste0('Spearman cor. with ',variable_list[[counter]][1])
          )
}
graphics.off()

