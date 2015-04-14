# This R script generates the quantile tables and the corresponding plots for the free parameters of different models of WCN
# Amir Shahmoradi, Sunday 9:57 PM, March 29, 2015, Wilke Lab, iCMB, UT Austin
# Note that R script input_data.r must be first sourced in order to use this script.

setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_best_definition/analysis/src')

dflist = list( wcneSC_bfSC   = exp_wcneSC_bfSC  
             , wcneSC_r4sJC  = exp_wcneSC_r4sJC 
             , wcneSC_seqent = exp_wcneSC_seqent
             , wcneSC_ddgent = exp_wcneSC_ddgent
             , wcngSC_bfSC   = exp_wcngSC_bfSC  
             , wcngSC_r4sJC  = exp_wcngSC_r4sJC 
             , wcngSC_seqent = exp_wcngSC_seqent
             , wcngSC_ddgent = exp_wcngSC_ddgent
             , wcnhSC_bfSC   = exp_wcnhSC_bfSC  
             , wcnhSC_r4sJC  = exp_wcnhSC_r4sJC 
             , wcnhSC_seqent = exp_wcnhSC_seqent
             , wcnhSC_ddgent = exp_wcnhSC_ddgent
             , wcnpSC_bfSC   = exp_wcnpSC_bfSC  
             , wcnpSC_r4sJC  = exp_wcnpSC_r4sJC 
             , wcnpSC_seqent = exp_wcnpSC_seqent
             , wcnpSC_ddgent = exp_wcnpSC_ddgent
             , wcneCA_bfSC   = exp_wcneCA_bfCA  
             , wcneCA_r4sJC  = exp_wcneCA_r4sJC 
             , wcneCA_seqent = exp_wcneCA_seqent
             , wcneCA_ddgent = exp_wcneCA_ddgent
             , wcngCA_bfSC   = exp_wcngCA_bfCA
             , wcngCA_r4sJC  = exp_wcngCA_r4sJC 
             , wcngCA_seqent = exp_wcngCA_seqent
             , wcngCA_ddgent = exp_wcngCA_ddgent
             , wcnhCA_bfSC   = exp_wcnhCA_bfCA  
             , wcnhCA_r4sJC  = exp_wcnhCA_r4sJC 
             , wcnhCA_seqent = exp_wcnhCA_seqent
             , wcnhCA_ddgent = exp_wcnhCA_ddgent
             , wcnpCA_bfSC   = exp_wcnpCA_bfCA  
             , wcnpCA_r4sJC  = exp_wcnpCA_r4sJC 
             , wcnpCA_seqent = exp_wcnpCA_seqent
             , wcnpCA_ddgent = exp_wcnpCA_ddgent
             )

counter = 0
for (dataframe in dflist)
{
  counter = counter + 1
  
  ###############################################################################
  ###############################################################################  
  # First write out the quantiles for the raw values of Spearman correlation strengths:
  quantiles = data.frame()
  for (i in names(dataframe[,-1])){
     row = data.frame(parameter     = dataframe[[i]][1],
                      mean_sp       = mean(dataframe[[i]][-1]),
                      median_sp     = quantile(as.vector(dataframe[[i]][-1]), probs = 0.50),
                      quantile05_sp = quantile(as.vector(dataframe[[i]][-1]), probs = 0.05),
                      quantile25_sp = quantile(as.vector(dataframe[[i]][-1]), probs = 0.25),
                      quantile75_sp = quantile(as.vector(dataframe[[i]][-1]), probs = 0.75),
                      quantile95_sp = quantile(as.vector(dataframe[[i]][-1]), probs = 0.95),
                      stdev_sp      = sd(as.vector(dataframe[[i]][-1]))
                      )
     quantiles = rbind(quantiles,row)
  }
  rownames(quantiles) = NULL
  write.csv( quantiles, file = paste0('../tables/get_quantiles/sp_raw/',names(dflist)[[counter]],'_quantiles.csv'), row.names=F)
  
  # Now generate the plot:
  model = substr(names(dflist)[[counter]],start=4,stop=4)
  atom  = substr(names(dflist)[[counter]],start=5,stop=6)
  var   = substr(names(dflist)[[counter]],start=8,stop=9)
  if (model == "e") {main_lab = "Exponential WCN"; xlab = 'Exponential Mean [Angstroms]'}
  if (model == "g") {main_lab = "Gaussian WCN"; xlab = 'Scale parameter [Angstroms]'}
  if (model == "h") {main_lab = "Heaviside WCN"; xlab = 'Cutoff Distance [Angstroms]'}
  if (model == "p") {main_lab = "Power-law WCN"; xlab = 'Power-law Exponent'; quantiles = quantiles[!quantiles$parameter==0.0,] }
  if (var == "bf") ylab = expression( paste( "Spearman ", rho ," : WCN - B factor" ) )
  if (var == "dd") ylab = expression( paste( "Spearman ", rho ," : WCN - ddG Rate" ) )
  if (var == "r4") ylab = expression( paste( "Spearman ", rho ," : WCN - r4sJC" ) )
  if (var == "se") ylab = expression( paste( "Spearman ", rho ," : WCN - Seq. Entropy" ) )
  
  #Plot Average raw Spearman correlation vs. value of free parameter
  pdf( paste0("../figures/get_quantiles/sp_raw/spcor_",names(dflist)[[counter]],".pdf"), width=4.5, height=4, useDingbats=FALSE )
  #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(quantiles$parameter,quantiles$mean_sp, type = 'l', ylim=c(-0.9,0.9), xlab=xlab, ylab=ylab, main=main_lab)
  polygon(c(quantiles$parameter,rev(quantiles$parameter)),c(quantiles$quantile25_sp,rev(quantiles$quantile75_sp)),col = "green", border = FALSE)
  lspl = smooth.spline(quantiles$parameter,quantiles$quantile25_sp)
  uspl = smooth.spline(quantiles$parameter,quantiles$quantile75_sp)
  lines(quantiles$parameter,quantiles$mean_sp,lwd=2)
  lines(quantiles$parameter,quantiles$median_sp,lwd=2,lty=2)
  #lines(predict(lspl,quantiles$parameter),lty=2,lwd=2,col='red')
  #lines(predict(uspl,quantiles$parameter),lty=2,lwd=2,col='red')
  graphics.off()


  ###############################################################################
  ###############################################################################
  # Now write out the quantiles for the absolute values of Spearman correlation strengths:
  quantiles = data.frame()
  for (i in names(dataframe[,-1])){
     row = data.frame(parameter     = dataframe[[i]][1],
                      mean_sp       = mean(abs(dataframe[[i]][-1])),
                      median_sp     = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.50),
                      quantile05_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.05),
                      quantile25_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.25),
                      quantile75_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.75),
                      quantile95_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.95),
                      stdev_sp      = sd(as.vector(dataframe[[i]][-1]))
                      )
     quantiles = rbind(quantiles,row)
  }
  rownames(quantiles) = NULL
  write.csv( quantiles, file = paste0('../tables/get_quantiles/sp_abs/',names(dflist)[[counter]],'_quantiles.csv'), row.names=F)

  # Now generate the plot:
  model = substr(names(dflist)[[counter]],start=4,stop=4)
  atom  = substr(names(dflist)[[counter]],start=5,stop=6)
  var   = substr(names(dflist)[[counter]],start=8,stop=9)
  if (model == "e") {main_lab = "Exponential WCN"; xlab = 'Exponential Mean [Angstroms]'}
  if (model == "g") {main_lab = "Gaussian WCN"; xlab = 'Scale parameter [Angstroms]'}
  if (model == "h") {main_lab = "Heaviside WCN"; xlab = 'Cutoff Distance [Angstroms]'}
  if (model == "p") {main_lab = "Power-law WCN"; xlab = 'Power-law Exponent'; quantiles = quantiles[!quantiles$parameter==0.0,] }
  if (var == "bf") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - B factor" ) )
  if (var == "dd") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - ddG Rate" ) )
  if (var == "r4") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - r4sJC" ) )
  if (var == "se") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - Seq. Entropy" ) )
  
  #Plot Average Absolute Spearman correlation vs. value of free parameter
  pdf( paste0("../figures/get_quantiles/sp_abs/spcor_",names(dflist)[[counter]],".pdf"), width=4.5, height=4, useDingbats=FALSE )
  #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(quantiles$parameter,quantiles$mean_sp, type = 'l', ylim=c(0.0,0.9), xlab=xlab, ylab=ylab, main=main_lab)
  polygon(c(quantiles$parameter,rev(quantiles$parameter)),c(quantiles$quantile25_sp,rev(quantiles$quantile75_sp)),col = "green", border = FALSE)
  lspl = smooth.spline(quantiles$parameter,quantiles$quantile25_sp)
  uspl = smooth.spline(quantiles$parameter,quantiles$quantile75_sp)
  lines(quantiles$parameter,quantiles$mean_sp,lwd=2)
  lines(quantiles$parameter,quantiles$median_sp,lwd=2,lty=2)
  lines(predict(lspl,quantiles$parameter),lty=2,lwd=2,col='red')
  lines(predict(uspl,quantiles$parameter),lty=2,lwd=2,col='red')
  graphics.off()


  ###############################################################################
  ###############################################################################
  #Plot Standard deviation of Spearman correlations vs. value of free parameters:
  if (var == "bf") ylab = expression( paste( "St. Dev. Spearman ", rho ," : WCN - B factor" ) )
  if (var == "dd") ylab = expression( paste( "St. Dev. Spearman ", rho ," : WCN - ddG Rate" ) )
  if (var == "r4") ylab = expression( paste( "St. Dev. Spearman ", rho ," : WCN - r4sJC" ) )
  if (var == "se") ylab = expression( paste( "St. Dev. Spearman ", rho ," : WCN - Seq. Entropy" ) )
  pdf( paste0("../figures/get_quantiles/stdev_",names(dflist)[[counter]],".pdf"), width=4.5, height=4, useDingbats=FALSE )
  #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(quantiles$parameter,quantiles$stdev_sp, type = 'l', lwd=2, xlab=xlab, ylab=ylab, main=main_lab)
  #lines(predict(lspl,quantiles$parameter),lty=2,lwd=2,col='red')
  #lines(predict(uspl,quantiles$parameter),lty=2,lwd=2,col='red')
  graphics.off()
}


##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
##################################################################################
# Make plots without main titles

counter = 0
for (dataframe in dflist)
{
  counter = counter + 1
  
  ###############################################################################
  ###############################################################################
  # write out the quantiles for the absolute values of Spearman correlation strengths:
  quantiles = data.frame()
  for (i in names(dataframe[,-1])){
     row = data.frame(parameter     = dataframe[[i]][1],
                      mean_sp       = mean(abs(dataframe[[i]][-1])),
                      median_sp     = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.50),
                      quantile05_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.05),
                      quantile25_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.25),
                      quantile75_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.75),
                      quantile95_sp = quantile(abs(as.vector(dataframe[[i]][-1])), probs = 0.95),
                      stdev_sp      = sd(as.vector(dataframe[[i]][-1]))
                      )
     quantiles = rbind(quantiles,row)
  }
  rownames(quantiles) = NULL
  
  # Now generate the plot:
  model = substr(names(dflist)[[counter]],start=4,stop=4)
  atom  = substr(names(dflist)[[counter]],start=5,stop=6)
  var   = substr(names(dflist)[[counter]],start=8,stop=9)
  if (model == "e") {xlab = 'Exponential Mean [ Angstroms ]'}
  if (model == "g") {xlab = 'Scale parameter [ Angstroms ]'}
  if (model == "h") {xlab = 'Cutoff Distance [ Angstroms ]'}
  if (model == "p") {xlab = 'Power-law Exponent'; quantiles = quantiles[!quantiles$parameter==0.0,] }
  if (var == "bf") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - B factor" ) )
  if (var == "dd") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - ddG Rate" ) )
  if (var == "r4") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - r4sJC" ) )
  if (var == "se") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - Seq. Entropy" ) )
  
  #Plot Average Absolute Spearman correlation vs. value of free parameter
  pdf( paste0("../figures/get_quantiles/screen_plots/spcor_",names(dflist)[[counter]],".pdf"), width=4.5, height=4, useDingbats=FALSE )
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot( quantiles$parameter,quantiles$mean_sp, type = 'l', ylim=c(0.0,0.9), xlab=xlab, ylab=ylab )
  polygon(c(quantiles$parameter,rev(quantiles$parameter)),c(quantiles$quantile25_sp,rev(quantiles$quantile75_sp)),col = "green", border = FALSE)
  lspl = smooth.spline(quantiles$parameter,quantiles$quantile25_sp)
  uspl = smooth.spline(quantiles$parameter,quantiles$quantile75_sp)
  lines(quantiles$parameter,quantiles$mean_sp,lwd=2)
  lines(quantiles$parameter,quantiles$median_sp,lwd=2,lty=2)
  lines(predict(lspl,quantiles$parameter),lty=2,lwd=2,col='red')
  lines(predict(uspl,quantiles$parameter),lty=2,lwd=2,col='red')
  graphics.off()
  
}