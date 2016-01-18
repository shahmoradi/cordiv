# Last updated by Amir Shahmoradi, Sunday 5:20 PM, January 18 2016, Wilke Lab, ICMB, UT Austin
# This code is exact same replicate of "best_structural_predictors_of_ER.r", except that all calculations are redone for the evolutionary rates calculated acording to JTT model.


setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
pdb_temp = data.frame()
best_structural_predictors_of_ER_JTT = data.frame()  # This data frame will contain correlations of select structural variables with r4sJC evolutionary rates, for each pdb file on each row
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JTT')]
  pdb_jec_ddg= res_prop_jec_ddg[res_prop_jec_ddg$pdb==substr(pdb,start=1,stop=4),] # c('rate.ddg.foldx')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  pdb_distance = res_prop_dfcSC[res_prop_dfcSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(pdb,seqent,ddgent))
                  , subset(pdb_jec, select = c(zr4s_JTT))
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
  r.rsa.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$rsa,method='sp')
  r.wcnSC.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$wcnSC,method='sp')
  r.wcnCA.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$wcnCA,method='sp')
  r.vareaSC.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$VSCarea,method='sp')
  r.vvolumeSC.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$VSCvolume,method='sp')
  #r.ddgent.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$ddgent,method='sp')
  r.ddgent.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$rate.ddg.foldx,method='sp')
  r.vareaCA.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$VCAarea,method='sp')
  r.vvolumeCA.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$VCAvolume,method='sp')
  r.bfSC.r4sJC = cor(pdb_temp$zr4s_JTT,pdb_temp$bfSC,method='sp')
  r.hbe.r4sJC = cor(pdb_temp$zr4s_JTT,abs(pdb_temp$hbe),method='sp')
  r.distance.r4sJC = cor(pdb_temp$zr4s_JTT,abs(pdb_temp$distance_normalized),method='sp')
  r.bfSC.wcnSC = cor(pdb_temp$wcnSC,pdb_temp$bfSC,method='sp')
  r.bfSC.distance = cor(pdb_temp$distance_normalized,pdb_temp$bfSC,method='sp')
  r.distance.wcnSC = cor(pdb_temp$distance_normalized,pdb_temp$wcnSC,method='sp')
  r.distance.vareaSC = cor(pdb_temp$distance_normalized,pdb_temp$VSCarea,method='sp')
  r.wcnSC.vareaSC = cor(pdb_temp$wcnSC,pdb_temp$VSCarea,method='sp')
  
  #row = data.frame( pdb = pdb, rsa = r.rsa.r4sJC, wcn = r.wcn.r4sJC, vareaSC = r.vareaSC.r4sJC, ddgent = r.ddgent.r4sJC )
  row = data.frame( pdb, r.wcnCA.r4sJC, r.wcnSC.r4sJC, r.rsa.r4sJC, r.vareaSC.r4sJC, r.vvolumeSC.r4sJC, r.vvolumeCA.r4sJC, r.ddgent.r4sJC , r.vareaCA.r4sJC , r.bfSC.r4sJC , r.hbe.r4sJC , r.distance.r4sJC, r.bfSC.distance, r.bfSC.wcnSC, r.distance.wcnSC , r.distance.vareaSC , r.wcnSC.vareaSC)
  best_structural_predictors_of_ER_JTT = rbind( best_structural_predictors_of_ER_JTT, row )
}

write.csv(best_structural_predictors_of_ER_JTT, file = "../tables/best_structural_predictors_of_ER_JTT.csv", row.names=F )

# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
hist.rsa = density(best_structural_predictors_of_ER_JTT$r.rsa.r4sJC)
hist.wcnSC = density(best_structural_predictors_of_ER_JTT$r.wcnSC.r4sJC)
hist.wcnCA = density(best_structural_predictors_of_ER_JTT$r.wcnCA.r4sJC)
hist.vareaSC = density(best_structural_predictors_of_ER_JTT$r.vareaSC.r4sJC)
hist.vvolumeSC = density(best_structural_predictors_of_ER_JTT$r.vvolumeSC.r4sJC)
hist.ddgent = density(best_structural_predictors_of_ER_JTT$r.ddgent.r4sJC)
hist.bfSC = density(best_structural_predictors_of_ER_JTT$r.bfSC.r4sJC)
hist.hbe = density(best_structural_predictors_of_ER_JTT$r.hbe.r4sJC)
hist.dist = density(best_structural_predictors_of_ER_JTT$r.distance.r4sJC)

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_structural_predictors_of_ER_JTT.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.vvolumeSC$x
    ,  hist.vvolumeSC$y
    ,   col = 'red'
    #,   xlim = c(-0.5,0.85)
    ,   xlim = c(-0.05,0.85)
    ,   ylim = c(0,5.5)
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
      , c(#"Distance from Geometrical Center of Protein",
        "Voronoi Cell Volume (SC)","1 / WCN (SC)", "1 / WCN (CA)", "B factor (SC)", "ddG Rate", "RSA")
      , col = c(#'grey',
        'red','black','black','cyan2','green','blue')
      , lty = c(#1,
        1,1,2,1,1,1)
      , lwd = 2
      , bty = 'n'
      , cex = 0.9
      )

graphics.off()


### The following plot is specifically generated for the manuscript:
# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_structural_predictors_of_ER_JTT_limited.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.vvolumeSC$x
       ,  hist.vvolumeSC$y
       ,   col = 'red'
       #,   xlim = c(-0.5,0.85)
       ,   xlim = c(0.1,0.9)
       ,   ylim = c(0,5.5)
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
lines( hist.rsa$x
       , hist.rsa$y
       , col = 'blue'
       , lwd  = 2
)
legend( 'topleft'
        , c("Voronoi Cell Volume (SC)","1 / WCN (SC)", "1 / WCN (CA)", "RSA")
        , col = c('red','black','black','blue')
        , lty = c(1,1,2,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)

graphics.off()


# Generate screen plots for the manuscript:
source("best_structural_predictors_of_ER_JTT_given_VSCvolume.r")

pdf( "../figures/best_structural_predictors_of_ER_JTT_screen.pdf", width=9, height=4, useDingbats=FALSE )
split.screen(c(1,2))

screen(1)

par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.vvolumeSC$x
       ,  hist.vvolumeSC$y
       ,   col = 'red'
       ,   xlim = c(0.1,0.9)
       ,   ylim = c(0,5.5)
       ,   type = 'l'
       ,   lwd  = 2 
       ,   xlab = "Correlation with ER"
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
lines( hist.rsa$x
       , hist.rsa$y
       , col = 'blue'
       , lwd  = 2
)
legend( 'topleft'
        , c("Voronoi Cell Volume (SC)","(-1) x WCN (SC)", "(-1) x WCN (CA)", "RSA")
        , col = c('red','black','black','blue')
        , lty = c(1,1,2,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)
mtext('A', side = 3, at=-0.04, font=2, cex=1.2)

screen(2)

par( mai=c(0.65, 0.65, 0.2, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  -5
       ,  -5
       ,   col = 'red'
       #,   xlim = c(-0.5,0.85)
       ,   xlim = c(-0.3,0.7)
       ,   ylim = c(0,4.)
       #,   col=colors[1]
       #,   ylim=c(0,7)
       #,   border = colors[1]
       #,   lty = 0
       ,   type = 'l'
       ,   lwd  = 2 
       #,   main = 'Correlations with Evolutionary Rates'
       #,   xlab = expression(paste('Spearman Cor. with Evolutionary Rates ',rho))
       ,   xlab = "Partial Correlation with ER given Cell Volume"
       ,   ylab = 'Relative Frequency'
)
lines( -hist.VwcnSC$x
       , hist.VwcnSC$y
       , col = 'black'
       , lwd = 2
)
lines( -hist.VwcnCA$x
       , hist.VwcnCA$y
       , col = 'black'
       , lwd = 2
       , lty = 2
)
lines( hist.Vrsa$x
       , hist.Vrsa$y
       , col = 'blue'
       , lwd  = 2
)

legend( 'topleft'
        , c("(-1) x WCN (SC)", "(-1) x WCN (CA)", "RSA")
        , col = c('black','black','blue')
        , lty = c(1,2,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)
mtext('B', side = 3, at=-0.475, font=2, cex=1.2)

close.screen(all = TRUE)

graphics.off()












# Now summarize the partial data over the entire dataset:

best_structural_predictors_of_ER_JTT = read.csv(file = "../tables/best_structural_predictors_of_ER_JTT.csv", header=TRUE )

data_subset = best_structural_predictors_of_ER_JTT[,-1]  # remove pdb column from the data set to loop over the columns below
best_structural_predictors_of_ER_JTT_summary = data.frame()
for (correlation in names(data_subset))
{
  row = data.frame( correlation = correlation
                    , mean = mean(best_structural_predictors_of_ER_JTT[[correlation]])
                    , median = median(best_structural_predictors_of_ER_JTT[[correlation]])
                    , min = min(best_structural_predictors_of_ER_JTT[[correlation]])
                    , quantile05 = quantile(best_structural_predictors_of_ER_JTT[[correlation]], probs=0.05)
                    , quantile25 = quantile(best_structural_predictors_of_ER_JTT[[correlation]], probs=0.25)
                    , quantile75 = quantile(best_structural_predictors_of_ER_JTT[[correlation]], probs=0.75)
                    , quantile95 = quantile(best_structural_predictors_of_ER_JTT[[correlation]], probs=0.95)
                    , max = max(best_structural_predictors_of_ER_JTT[[correlation]])
                    , sd = sd(best_structural_predictors_of_ER_JTT[[correlation]])
  )
  best_structural_predictors_of_ER_JTT_summary = rbind(best_structural_predictors_of_ER_JTT_summary,row)  
}
write.csv(best_structural_predictors_of_ER_JTT_summary, file = "../tables/best_structural_predictors_of_ER_JTT_summary.csv", row.names=F )






max(best_structural_predictors_of_ER_JTT$r.wcnSC.r4sJC)-min(best_structural_predictors_of_ER_JTT$r.wcnSC.r4sJC)
max(best_structural_predictors_of_ER_JTT$r.vareaSC.r4sJC)-min(best_structural_predictors_of_ER_JTT$r.vareaSC.r4sJC)
sd(best_structural_predictors_of_ER_JTT$r.vareaSC.r4sJC)
sd(best_structural_predictors_of_ER_JTT$r.wcnSC.r4sJC)

# Now do some paired t-tests:

temp1 = data.frame(var = 'wcnSC', r = best_structural_predictors_of_ER_JTT$r.wcnSC.r4sJC)
temp2 = data.frame(var = 'vareaSC', r = best_structural_predictors_of_ER_JTT$r.vareaSC.r4sJC)
temp = rbind(temp1,temp2)
temp$var = factor(temp$var)
ttest = pairwise.t.test(temp$r,temp$var)

install.packages('reshape')
library('reshape')
best_structural_predictors_of_ER_JTT_long = reshape(best_structural_predictors_of_ER_JTT, direction='long', varying=colnames(best_structural_predictors_of_ER_JTT)[2:ncol(best_structural_predictors_of_ER_JTT)], idvar=c('pdb'), v.names='value', timevar='variable', times=colnames(best_structural_predictors_of_ER_JTT)[2:ncol(best_structural_predictors_of_ER_JTT)])
best_structural_predictors_of_ER_JTT_long$variable = factor(best_structural_predictors_of_ER_JTT_long$variable)
best_structural_predictors_of_ER_JTT_long$value[best_structural_predictors_of_ER_JTT_long$variable=='r.wcnCA.r4sJC'] = -best_structural_predictors_of_ER_JTT_long$value[best_structural_predictors_of_ER_JTT_long$variable=='r.wcnCA.r4sJC']
best_structural_predictors_of_ER_JTT_long$value[best_structural_predictors_of_ER_JTT_long$variable=='r.wcnSC.r4sJC'] = -best_structural_predictors_of_ER_JTT_long$value[best_structural_predictors_of_ER_JTT_long$variable=='r.wcnSC.r4sJC']

ttest = pairwise.t.test(best_structural_predictors_of_ER_JTT_long$value
                       ,best_structural_predictors_of_ER_JTT_long$variable
                       ,p.adjust.method = "hochberg"
                       ,paired=TRUE
                       )

pvalues = as.data.frame(ttest$p.value)

write.csv(pvalues, file='../tables/pairwise_t_test.csv')

