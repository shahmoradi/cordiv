# This R script generates the quantile tables and the corresponding plots for the free parameters of different models of WCN
# Amir Shahmoradi, Sunday 9:57 PM, March 29, 2015, Wilke Lab, iCMB, UT Austin
# Note that R script input_data.r must be first sourced in order to use this script.

setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_best_definition/analysis/src')

dflist = list( wcneSC_bfSC   = exp_wcneSC_bfSC
             , wcngSC_bfSC   = exp_wcngSC_bfSC
             , wcnhSC_bfSC   = exp_wcnhSC_bfSC
             , wcnpSC_bfSC   = exp_wcnpSC_bfSC
             , wcneCA_bfCA   = exp_wcneCA_bfCA
             , wcngCA_bfCA   = exp_wcngCA_bfCA
             , wcnhCA_bfCA   = exp_wcnhCA_bfCA
             , wcnpCA_bfCA   = exp_wcnpCA_bfCA
             , wcneSC_r4sJC  = exp_wcneSC_r4sJC
             , wcngSC_r4sJC  = exp_wcngSC_r4sJC
             , wcnhSC_r4sJC  = exp_wcnhSC_r4sJC
             , wcnpSC_r4sJC  = exp_wcnpSC_r4sJC
             , wcneCA_r4sJC  = exp_wcneCA_r4sJC
             , wcngCA_r4sJC  = exp_wcngCA_r4sJC
             , wcnhCA_r4sJC  = exp_wcnhCA_r4sJC
             , wcnpCA_r4sJC  = exp_wcnpCA_r4sJC
             , wcneSC_seqent = exp_wcneSC_seqent
             , wcngSC_seqent = exp_wcngSC_seqent
             , wcnhSC_seqent = exp_wcnhSC_seqent
             , wcnpSC_seqent = exp_wcnpSC_seqent
             , wcneCA_seqent = exp_wcneCA_seqent
             , wcngCA_seqent = exp_wcngCA_seqent
             , wcnhCA_seqent = exp_wcnhCA_seqent
             , wcnpCA_seqent = exp_wcnpCA_seqent
             , wcneSC_ddgent = exp_wcneSC_ddgent
             , wcngSC_ddgent = exp_wcngSC_ddgent
             , wcnhSC_ddgent = exp_wcnhSC_ddgent
             , wcnpSC_ddgent = exp_wcnpSC_ddgent
             , wcneCA_ddgent = exp_wcneCA_ddgent
             , wcngCA_ddgent = exp_wcngCA_ddgent
             , wcnhCA_ddgent = exp_wcnhCA_ddgent
             , wcnpCA_ddgent = exp_wcnpCA_ddgent
             )
counter = 0
for (dataframe in dflist){
  counter = counter + 1
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
  write.csv( quantiles, file = paste0('../tables/',names(dflist)[[counter]],'_quantiles.csv'), row.names=F)
  
  # Now generate the plot:
  model = substr(names(dflist)[[counter]],start=4,stop=4)
  atom  = substr(names(dflist)[[counter]],start=5,stop=6)
  var   = substr(names(dflist)[[counter]],start=8,stop=9)
  if (model == "e") {main_lab = "Exponential WCN"; xlab = 'Exponential Mean [Angstroms]'}
  if (model == "g") {main_lab = "Gaussian WCN"; xlab = 'Standard Deviation [Angstroms]'}
  if (model == "h") {main_lab = "Heaviside WCN"; xlab = 'Cutoff Distance [Angstroms]'}
  if (model == "p") {main_lab = "Power-law WCN"; xlab = 'Power-law Exponent'; quantiles = quantiles[!quantiles$parameter==0.0,] }
  if (var == "bf") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - B factor" ) )
  if (var == "dd") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - ddG Rate" ) )
  if (var == "r4") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - r4sJC" ) )
  if (var == "se") ylab = expression( paste( "Absolute Spearman ", rho ," : WCN - Seq. Entropy" ) )
  
  pdf( paste0("../figures/spcor_",names(dflist)[[counter]],".pdf"), width=4.5, height=4, useDingbats=FALSE )
  #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(quantiles$parameter,quantiles$mean_sp, type = 'l', ylim=c(0.,0.9), xlab=xlab, ylab=ylab, main=main_lab)
  polygon(c(quantiles$parameter,rev(quantiles$parameter)),c(quantiles$quantile25_sp,rev(quantiles$quantile75_sp)),col = "green", border = FALSE)
  lspl = smooth.spline(quantiles$parameter,quantiles$quantile25_sp)
  uspl = smooth.spline(quantiles$parameter,quantiles$quantile75_sp)
  lines(quantiles$parameter,quantiles$mean_sp,lwd=2)
  lines(quantiles$parameter,quantiles$median_sp,lwd=2,lty=2)
  lines(predict(lspl,quantiles$parameter),lty=2,lwd=2,col='red')
  lines(predict(uspl,quantiles$parameter),lty=2,lwd=2,col='red')
  graphics.off()
}
