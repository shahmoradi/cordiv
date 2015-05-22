# This R script is very similar to best_structural_predictors_of_ER.r , with the exception that here the PARTIAL correlations are measured while controlling for Voronoi cell volume.

# Last updated by Amir Shahmoradi, Thursday 4:25 PM, May 21 2015, Wilke Lab, ICMB, UT Austin

library(ppcor)
setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
pdb_temp = data.frame()
best_structural_predictors_of_ER_given_VSCvolume = data.frame()  # This data frame will contain correlations of select structural variables with r4sJC evolutionary rates, for each pdb file on each row
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
  pdb_jec_ddg= res_prop_jec_ddg[res_prop_jec_ddg$pdb==substr(pdb,start=1,stop=4),] # c('rate.ddg.foldx')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  pdb_distance = res_prop_dfcSC[res_prop_dfcSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(pdb,seqent,ddgent))
                  , subset(pdb_jec, select = c(zr4s_JC))
                  , subset(pdb_jec_ddg, select = c(rate.ddg.foldx))
                  , subset(pdb_hps, select = c(hpshh))
                  , subset(pdb_dssp, select = c(rsa,hbe))
                  , subset(pdb_wcn_bf, select = c(wcnSC,wcnCA,bfSC))
                  , subset(pdb_voroSC, select = c(VSCarea))
                  , subset(pdb_voroCA, select = c(VCAarea))
                  , subset(pdb_voroSC, select = c(VSCvolume))
                  , subset(pdb_voroCA, select = c(VCAvolume))
                  , subset(pdb_distance, select = c(distance_normalized))
                  )
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$rsa,pdb_temp$VSCvolume,method='sp'); r.rsa.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$wcnSC,pdb_temp$VSCvolume,method='sp'); r.wcnSC.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$wcnCA,pdb_temp$VSCvolume,method='sp'); r.wcnCA.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$VSCarea,pdb_temp$VSCvolume,method='sp'); r.vareaSC.r4sJC = x$estimate
 #x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$ddgent,pdb_temp$VSCvolume,method='sp'); r.ddgent.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$rate.ddg.foldx,pdb_temp$VSCvolume,method='sp'); r.ddgent.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$VCAarea,pdb_temp$VSCvolume,method='sp'); r.vareaCA.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$VCAvolume,pdb_temp$VSCvolume,method='sp'); r.vvolumeCA.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,pdb_temp$bfSC,pdb_temp$VSCvolume,method='sp'); r.bfSC.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,abs(pdb_temp$hbe),pdb_temp$VSCvolume,method='sp'); r.hbe.r4sJC = x$estimate
  x = pcor.test(pdb_temp$zr4s_JC,abs(pdb_temp$distance_normalized),pdb_temp$VSCvolume,method='sp'); r.distance.r4sJC = x$estimate
  x = pcor.test(pdb_temp$wcnSC,pdb_temp$bfSC,pdb_temp$VSCvolume,method='sp'); r.bfSC.wcnSC = x$estimate
  x = pcor.test(pdb_temp$distance_normalized,pdb_temp$bfSC,pdb_temp$VSCvolume,method='sp'); r.bfSC.distance = x$estimate
  x = pcor.test(pdb_temp$distance_normalized,pdb_temp$wcnSC,pdb_temp$VSCvolume,method='sp'); r.distance.wcnSC = x$estimate
  x = pcor.test(pdb_temp$distance_normalized,pdb_temp$VSCarea,pdb_temp$VSCvolume,method='sp'); r.distance.vareaSC = x$estimate
  x = pcor.test(pdb_temp$wcnSC,pdb_temp$VSCarea,pdb_temp$VSCvolume,method='sp'); r.wcnSC.vareaSC = x$estimate
  
  #row = data.frame( pdb = pdb, rsa = r.rsa.r4sJC, wcn = r.wcn.r4sJC, vareaSC = r.vareaSC.r4sJC, ddgent = r.ddgent.r4sJC )
  row = data.frame( pdb, r.wcnCA.r4sJC, r.wcnSC.r4sJC, r.rsa.r4sJC, r.vareaSC.r4sJC, r.vvolumeCA.r4sJC, r.ddgent.r4sJC , r.vareaCA.r4sJC , r.bfSC.r4sJC , r.hbe.r4sJC , r.distance.r4sJC, r.bfSC.distance, r.bfSC.wcnSC, r.distance.wcnSC , r.distance.vareaSC , r.wcnSC.vareaSC)
  best_structural_predictors_of_ER_given_VSCvolume = rbind( best_structural_predictors_of_ER_given_VSCvolume, row )
}

write.csv(best_structural_predictors_of_ER_given_VSCvolume, file = "../tables/best_structural_predictors_of_ER_given_VSCvolume.csv", row.names=F )

# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
hist.rsa = density(best_structural_predictors_of_ER_given_VSCvolume$r.rsa.r4sJC)
hist.wcnSC = density(best_structural_predictors_of_ER_given_VSCvolume$r.wcnSC.r4sJC)
hist.wcnCA = density(best_structural_predictors_of_ER_given_VSCvolume$r.wcnCA.r4sJC)
hist.vareaSC = density(best_structural_predictors_of_ER_given_VSCvolume$r.vareaSC.r4sJC)
hist.ddgent = density(best_structural_predictors_of_ER_given_VSCvolume$r.ddgent.r4sJC)
hist.bfSC = density(best_structural_predictors_of_ER_given_VSCvolume$r.bfSC.r4sJC)
hist.hbe = density(best_structural_predictors_of_ER_given_VSCvolume$r.hbe.r4sJC)
hist.dist = density(best_structural_predictors_of_ER_given_VSCvolume$r.distance.r4sJC)

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_structural_predictors_of_ER_given_VSCvolume.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  -5
    ,  -5
    ,   col = 'red'
    #,   xlim = c(-0.5,0.85)
    ,   xlim = c(-0.3,0.7)
    ,   ylim = c(0,5.)
    #,   col=colors[1]
    #,   ylim=c(0,7)
    #,   border = colors[1]
    #,   lty = 0
    ,   type = 'l'
    ,   lwd  = 2 
    #,   main = 'Correlations with Evolutionary Rates'
    #,   xlab = expression(paste('Spearman Cor. with Evolutionary Rates ',rho))
    ,   xlab = expression( paste('Correlation with Evolutionary Rates: Spearman ', rho ) )
    ,   ylab = 'Relative Frequency'
    )
lines( -hist.wcnSC$x
     , hist.wcnSC$y
     , col = 'black'
     , lwd = 2
     )
lines( -hist.wcnCA$x
       , hist.wcnCA$y
       , col = 'black'
       , lwd = 2
       , lty = 2
     )
lines( hist.bfSC$x
       , hist.bfSC$y
       , col = 'cyan2'
       , lwd = 2
     )
lines( hist.ddgent$x
       , hist.ddgent$y
       , col = 'green'
       , lwd  = 2
     )
lines( hist.rsa$x
       , hist.rsa$y
       , col = 'blue'
       , lwd  = 2
     )
#lines( hist.dist$x
#       , hist.dist$y
#       , col = 'grey'
#       , lwd  = 2
#     )

#lines( -hist.hbe$x
#     , hist.hbe$y
#     , col = 'grey'
#     , lwd = 2
#     )

#legend( 'topleft'
#      , c("Voronoi Cell Area","H-bond strength", "WCN (SC)", "WCN (CA)", "ddG Rate", "Bfactor", "RSA")
#      , col = c('red','grey','black','black','green','cyan2','blue')
#lines( -hist.hbe$x
#     , hist.hbe$y
#     , col = 'grey'
#     , lwd = 2
#     )

legend( 'topleft'
      , c(#"Distance from Geometrical Center of Protein","Voronoi Cell Volume (SC)",
        "1 / WCN (SC)", "1 / WCN (CA)", "B factor (SC)", "ddG Rate", "RSA")
      , col = c(#'grey','red',
        'black','black','cyan2','green','blue')
      , lty = c(#1,1,
        1,2,1,1,1)
      , lwd = 2
      , bty = 'n'
      , cex = 0.9
      )

graphics.off()





max(best_structural_predictors_of_ER_given_VSCvolume$r.wcnSC.r4sJC)-min(best_structural_predictors_of_ER_given_VSCvolume$r.wcnSC.r4sJC)
max(best_structural_predictors_of_ER_given_VSCvolume$r.vareaSC.r4sJC)-min(best_structural_predictors_of_ER_given_VSCvolume$r.vareaSC.r4sJC)
sd(best_structural_predictors_of_ER_given_VSCvolume$r.vareaSC.r4sJC)
sd(best_structural_predictors_of_ER_given_VSCvolume$r.wcnSC.r4sJC)

# Now do some paired t-tests:

#temp1 = data.frame(var = 'wcnSC', r = best_structural_predictors_of_ER_given_VSCvolume$r.wcnSC.r4sJC)
#temp2 = data.frame(var = 'vareaSC', r = best_structural_predictors_of_ER_given_VSCvolume$r.vareaSC.r4sJC)
#temp = rbind(temp1,temp2)
#temp$var = factor(temp$var)
#ttest = pairwise.t.test(temp$r,temp$var)

#install.packages('reshape')
#library('reshape')
#best_structural_predictors_of_ER_given_VSCvolume_long = reshape(best_structural_predictors_of_ER_given_VSCvolume, direction='long', varying=colnames(best_structural_predictors_of_ER_given_VSCvolume)[2:ncol(best_structural_predictors_of_ER_given_VSCvolume)], idvar=c('pdb'), v.names='value', timevar='variable', times=colnames(best_structural_predictors_of_ER_given_VSCvolume)[2:ncol(best_structural_predictors_of_ER_given_VSCvolume)])
#best_structural_predictors_of_ER_given_VSCvolume_long$variable = factor(best_structural_predictors_of_ER_given_VSCvolume_long$variable)

#ttest = pairwise.t.test(best_structural_predictors_of_ER_given_VSCvolume_long$value,best_structural_predictors_of_ER_given_VSCvolume_long$variable)

#pvalues = as.data.frame(ttest$p.value)

#write.csv(pvalues, file='../tables/pairwise_t_test_given_VSCvolume.csv')
