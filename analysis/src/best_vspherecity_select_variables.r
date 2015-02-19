
# This code is aimed at finding which definition of Voronoi cell volume (vvol) and area (varea) can best represent a residue. This is done by comparing the correlations of different vorvols with other residue variables.
# "select_variables" in the name of this code refers to the fact that I am comparing the performance of different vorvols to only a select number of important residue properties, such as ASA, RSA, seqent, ddGent, and I am ignoring the rest which are generally definition dependent, such as different definitions of B factors and varea based on the set of atom coordinates used.

# Last updated by Amir Shahmoradi, Friday 4:04 PM, February 13 2015, Wilke Lab, ICMB, UT Austin

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

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

varea_scors_all_pdbs = data.frame()
# The above dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
# Here scor stands for Sperman CORrelation.
varea_list = c('vareaSC','vareaAA','vareaCA','vareaCB','vareaN','vareaC','vareaO') #,'wcnCA','wcnN','wcnC','wcnO')
crd_list = c('SC','AA','CB','CA','N','C','O') #,'wcnCA','wcnN','wcnC','wcnO') # list of representative coordinates used to generate Voronoi cells.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('r4sJC')]
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,]  # c('hpskd','hpsww','hpshh')]
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  #pdb_bf    = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb]
  pdb_varea  = cbind( res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
                    , res_prop_voroAA[res_prop_voroAA$pdb==pdb, ]
                    , res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
                    , res_prop_voroCB[res_prop_voroCB$pdb==pdb, ]
                    , res_prop_voroN[res_prop_voroN$pdb==pdb, ]
                    , res_prop_voroC[res_prop_voroC$pdb==pdb, ]
                    , res_prop_voroO[res_prop_voroO$pdb==pdb, ] )
  
  pdb_varea  = data.frame( vareaSC = pdb_varea$VSCarea 
                         , vareaAA = pdb_varea$VAAarea 
                         , vareaCA = pdb_varea$VCAarea
                         , vareaCB = pdb_varea$VCBarea
                         , vareaN = pdb_varea$VNarea
                         , vareaC = pdb_varea$VCarea
                         , vareaO = pdb_varea$VOarea )

  pdb_temp = cbind( subset(pdb_jec, select = c(r4s_JC))
                  , subset(pdb_elj, select = c(seqent,ddgent))
                  , subset(pdb_hps, select  = c(hpshh))
                  # subset(pdb_dssp, select = c(asa,rsa,hbe))
                  , subset(pdb_dssp, select = c(rsa,hbe))
                  # subset(pdb_voroAA, select = c(resvol))
                  # subset(pdb_varea, select = -c(pdb,resnam,resnum,bfSC,bfAA,bfN,bfCA,bfC,bfO,bfCB)),
                  # subset(pdb_voroAA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,VAAnvertices,VAAnedges,VAAvolume_change)),
                  # subset(pdb_voroCA, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VCAnvertices,VCAnedges,VCAvolume_change)),
                  # subset(pdb_voroSC, select = -c(pdb,resnam,resnum,sizeSC,sizeAA,resvol,VSCnvertices,VSCnedges,VSCvolume_change))
                  )
  
  pdb_varea_long = reshape(pdb_varea, ids = rownames(pdb_varea), varying = colnames(pdb_varea), v.names = 'value', timevar = 'variable', times = colnames(pdb_varea), direction = 'long')
  pdb_varea_long$variable = factor(pdb_varea_long$variable)
  pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
  pdb_long$variable = factor(pdb_long$variable)
  
  for (varea in levels(pdb_varea_long$variable))
  {
    varea_data = pdb_varea_long[pdb_varea_long$variable == varea,]
    for (variable in levels(pdb_long$variable))
    {
      variable_data = pdb_long[pdb_long$variable == variable,]
      x = cor.test( varea_data$value, variable_data$value, method='spearman', na.action="na.omit" )
      r = x$estimate
      p = x$p.value
      
      row = data.frame(pdb = pdb, varea = varea, variable = variable, value = r)
      varea_scors_all_pdbs = rbind(varea_scors_all_pdbs,row)
    }
  }
}
write.csv( varea_scors_all_pdbs, "../tables/best_varea/select_variables/varea_scors_all_pdbs.csv", row.names=F )



# Now summarize all Spearman correlations over the entire dataset

varea_scors_all_pdbs = read.csv( "../tables/best_varea/select_variables/varea_scors_all_pdbs.csv", header = TRUE )
varea_scors_all_pdbs$variable = factor(varea_scors_all_pdbs$variable)
varea_scors_all_pdbs$varea = factor(varea_scors_all_pdbs$varea)


for (varea in levels(varea_scors_all_pdbs$varea))
{
  temp_data_varea = varea_scors_all_pdbs[varea_scors_all_pdbs$varea == varea,]
  temp_data_varea$variable = factor(temp_data_varea$variable)
  
  varea_scors_summary = data.frame()
  for (variable in levels(temp_data_varea$variable))
  {
    temp_data = temp_data_varea[temp_data_varea$variable == variable,]
    if (length(temp_data$pdb) != npdbs)
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
    varea_scors_summary = rbind(varea_scors_summary,row)
  }
  row.names(varea_scors_summary) = c()
  filename = paste0("../tables/best_varea/select_variables/",varea,'_scors_summary.csv')
  write.csv(varea_scors_summary, filename, row.names=F )
  cat (varea, filename, '\n')
}

# Now compare the three varea scors with each other:

for (i in 1:(length(varea_list)-1))
{
  filename = paste0("../tables/best_varea/select_variables/",varea_list[[i]][1],'_scors_summary.csv')
  varea_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(varea_list))
  {
    filename = paste0("../tables/best_varea/select_variables/",varea_list[[j]][1],'_scors_summary.csv')
    varea_scors_summary2 = read.csv( filename, header = T )
    difference = data.frame(variable      = varea_scors_summary1$variable,
                            mean_diff     = abs(varea_scors_summary1$mean)   - abs(varea_scors_summary2$mean),
                            median_diff   = abs(varea_scors_summary1$median) - abs(varea_scors_summary2$median),
                            sd_diff       = abs(varea_scors_summary1$sd)     - abs(varea_scors_summary2$sd),
                            min_diff      = abs(varea_scors_summary1$min)    - abs(varea_scors_summary2$min),
                            q1_diff       = abs(varea_scors_summary1$q1)     - abs(varea_scors_summary2$q1),
                            q3_diff       = abs(varea_scors_summary1$q3)     - abs(varea_scors_summary2$q3),
                            max_diff      = abs(varea_scors_summary1$max)    - abs(varea_scors_summary2$max)
                            )
    row.names(difference) = c()
    filename = paste0('../tables/best_varea/select_variables/diff_',varea_list[[i]][1],'_',varea_list[[j]][1],'.csv')
    cat (filename, '\n')
    write.csv( difference, filename, row.names=F )
  }
}

# Make plots of ABS(varea) vs. ABS(varea) for comparison:

x = -2:2
colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2') #, 'darkred', 'darkgreen', 'bisque2')
#labels = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')
labels = c('ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Evol. Rates', 'Residue Volume', 'RSA', 'Seq. Entropy')

for (i in 1:(length(varea_list)-1))
{
  filename = paste0("../tables/best_varea/select_variables/",varea_list[[i]][1],'_scors_summary.csv')
  varea_scors_summary1 = read.csv( filename, header = T )
  for (j in (i+1):length(varea_list))
  {
    filename = paste0("../tables/best_varea/select_variables/",varea_list[[j]][1],'_scors_summary.csv')
    varea_scors_summary2 = read.csv( filename, header = T )
    filename = paste0('../figures/best_varea/select_variables/abs(',varea_list[[i]][1],'_',varea_list[[j]][1],').pdf')
    cat (filename, '\n')
    pdf( filename, width=4.5, height=4, useDingbats=FALSE )
    par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    plot(-999,
         #xaxt='n',yaxt='n',bty='n',pch='',
    #plot(abs(varea_scors_summary1$median),
    #     abs(varea_scors_summary2$median),
         xlab = paste0( 'absolute median cor. with ',varea_list[[i]][1] ),
         ylab = paste0( 'absolute median cor. with ',varea_list[[j]][1] ),
         xlim = c(0,1.0),
         ylim = c(0,1.0)
         )
    points( abs(varea_scors_summary1$median),
            abs(varea_scors_summary2$median), pch=19, col=colors[1:length(varea_scors_summary2$median)])
    lines( x, x, col='red' )
    legend( -0.02, 1.05, labels[1:7], pch=19, col=colors[1:7], bty='n', cex=0.9)
    #legend( 0.7, 0.1, labels[4:6], pch=19, col=colors[4:6], bty='n', cex=0.9)
    graphics.off()
  }
}

# Now create box plots

varea_scors_all_pdbs = read.csv( "../tables/best_varea/select_variables/varea_scors_all_pdbs.csv", header = TRUE )
varea_scors_all_pdbs$variable = factor(varea_scors_all_pdbs$variable)
varea_scors_all_pdbs$varea = factor(varea_scors_all_pdbs$varea)

#variable_list = c('ASA', 'ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Residue Volume', 'RSA', 'Seq. Entropy')
variable_list = c('ddG Entropy', 'H-bond energy', 'Hydrophobicity', 'Evol. Rates', 'Residue Volume', 'RSA', 'Seq. Entropy')
counter = 0
for (variable in levels(varea_scors_all_pdbs$variable))
{ 
  counter = counter + 1
  temp_data_variable_long = varea_scors_all_pdbs[varea_scors_all_pdbs$variable == variable,c('pdb','varea','value')]
  temp_data_variable = reshape(temp_data_variable_long, timevar = 'varea', idvar = 'pdb', direction = 'wide')
  temp_data_variable = subset (temp_data_variable, select = -c(pdb))
  #colnames(temp_data_variable) = levels(varea_scors_all_pdbs$varea)
  colnames(temp_data_variable) = substr(levels(varea_scors_all_pdbs$varea),start=6,stop=7)  # remove 'varea' from the names of the variables. It is unnecessary and ugly in the figures. Only the the atom names will remain.
  #temp_data_variable = temp_data_variable[varea_list]
  temp_data_variable = temp_data_variable[crd_list]
  filename = paste0('../figures/best_varea/select_variables/boxplot_',variable,'.pdf')
  pdf( filename, width=6, height=4, useDingbats=FALSE )
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  boxplot(temp_data_variable,
          xlab = 'representative Voronoi Cell Area (varea)',
          ylab = paste0('Spearman cor. with ',variable_list[[counter]][1])
  )
  graphics.off()
}

# Now create box plots all in one figure

varea_scors_all_pdbs = read.csv( "../tables/best_varea/select_variables/varea_scors_all_pdbs.csv", header = TRUE )
varea_scors_all_pdbs$variable = factor(varea_scors_all_pdbs$variable)
varea_scors_all_pdbs$varea = factor(varea_scors_all_pdbs$varea)

#variable_list = c('Seq. Entropy','ddG Entropy','RSA','ASA','Hydrophobicity','H-bond Energy')
#variable_names = c('seqent','ddgent','rsa','asa','hpshh','hbe')
variable_list = c('Evol. Rates','Seq. Entropy','ddG Entropy','RSA','Hydrophobicity','H-bond Energy')
variable_names = c('r4s_JC','seqent','ddgent','rsa','hpshh','hbe')
#varea_list = c('vareaSC','vareaAA','vareaCA') #,'wcnCA','wcnN','wcnC','wcnO')
counter = 0
filename = paste0('../figures/best_varea/select_variables/boxplot_varea_all_in_one.pdf')
pdf( filename, width=15, height=12, useDingbats=FALSE )
split.screen(c(3,2))
for (variable in variable_names)
{
  counter = counter + 1
  screen(counter)
  temp_data_variable_long = varea_scors_all_pdbs[varea_scors_all_pdbs$variable == variable,c('pdb','varea','value')]
  #temp_data_variable_long$varea = substr(temp_data_variable_long$varea,start=6,stop=7)  # remove 'varea' from the names of the variables. It is unnecessary and ugly in the figures. Only the the atom names will remain.
  temp_data_variable = reshape(temp_data_variable_long, timevar = 'varea', idvar = 'pdb', direction = 'wide')
  temp_data_variable = subset (temp_data_variable, select = -c(pdb))
  colnames(temp_data_variable) = substr(levels(varea_scors_all_pdbs$varea),start=6,stop=7)  # remove 'varea' from the names of the variables. It is unnecessary and ugly in the figures. Only the the atom names will remain.
  temp_data_variable = temp_data_variable[crd_list]
  #temp_data_variable = temp_data_variable[varea_list]
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  boxplot(temp_data_variable,
          xlab = 'representative Voronoi Cell Area (varea)',
          ylab = paste0('Spearman cor. with ',variable_list[[counter]][1]),
          cex.axis = 1.3,
          cex.lab = 1.3
          #, ylim = c(0.05,-1.)
          )
}
graphics.off()