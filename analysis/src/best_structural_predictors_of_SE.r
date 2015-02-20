# This code is aimed at finding which structural property best predicts the Sequence Entropy.

# Last updated by Amir Shahmoradi, Thursday 5:48 PM, February 19 2015, Wilke Lab, ICMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
pdb_temp = data.frame()
SEscors = data.frame()  # This data frame will contain correlations of select structural variables with Sequence Entropy, for each pdb file on eaxch row
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voro   = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(pdb,seqent,ddgent))
                  , subset(pdb_jec, select = c(zr4s_JC))
                  , subset(pdb_hps, select = c(hpshh))
                  , subset(pdb_dssp, select = c(rsa,hbe))
                  , subset(pdb_wcn_bf, select = c(wcnSC,wcnCA,bfSC))
                  , subset(pdb_voro, select = c(VSCarea))
                  )
  r.rsa.seqent = cor(pdb_temp$seqent,pdb_temp$rsa,method='sp')
  r.wcnSC.seqent = cor(pdb_temp$seqent,pdb_temp$wcnSC,method='sp')
  r.wcnCA.seqent = cor(pdb_temp$seqent,pdb_temp$wcnCA,method='sp')
  r.varea.seqent = cor(pdb_temp$seqent,pdb_temp$VSCarea,method='sp')
  r.ddgent.seqent = cor(pdb_temp$seqent,pdb_temp$ddgent,method='sp')
  
  #row = data.frame( pdb = pdb, rsa = r.rsa.seqent, wcn = r.wcn.seqent, varea = r.varea.seqent, ddgent = r.ddgent.seqent )
  row = data.frame( pdb, r.wcnSC.seqent, r.wcnCA.seqent, r.rsa.seqent, r.varea.seqent, r.ddgent.seqent )
  SEscors = rbind( SEscors, row )
}


# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
hist.rsa = density(SEscors$r.rsa.seqent)
hist.wcnSC = density(SEscors$r.wcnSC.seqent)
hist.wcnCA = density(SEscors$r.wcnCA.seqent)
hist.varea = density(SEscors$r.varea.seqent)
hist.ddgent = density(SEscors$r.ddgent.seqent)

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_structural_predictors_of_SE.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.rsa$x
    ,  hist.rsa$y
    ,   col = 'green'
    ,   xlim = c(0.0,0.8)
    ,   ylim = c(0,5.5)
    #,   col=colors[1]
    #,   ylim=c(0,7)
    #,   border = colors[1]
    #,   lty = 0
    ,   type = 'l'
    ,   lwd  = 2 
    #,   main = 'Correlations with Evolutionary Rates'
    #,   xlab = expression(paste('Absolute Spearman Cor. with Evolutionary Rates ',rho))
    ,   xlab = 'Absolute Spearman Cor. with Sequence Entropy'
    ,   ylab = 'frequency'
    )
lines( abs(hist.ddgent$x)
     , abs(hist.ddgent$y)
     , col = 'blue'
     , lwd  = 2
     )
lines( abs(hist.varea$x)
     , abs(hist.varea$y)
     , col = 'red'
     , lwd  = 2
     )
lines( abs(hist.wcnSC$x)
     , abs(hist.wcnSC$y)
     , col = 'black'
     , lwd = 2
     )
lines( abs(hist.wcnCA$x)
     , abs(hist.wcnCA$y)
     , col = 'black'
     , lwd = 2
     , lty = 2
     )

legend( 'topleft'
      , c("Voronoi Cell Area", "WCN (SC)", "WCN (CA)", "ddG Rate", "RSA")
      , col = c('red','black','black','green','blue')
      , lty = c(1,1,2,1,1)
      , lwd = 2
      , bty = 'n'
      , cex = 0.9
      )

graphics.off()
