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
  pdb_jec_ddg= res_prop_jec_ddg[res_prop_jec_ddg$pdb==substr(pdb,start=1,stop=4),] # c('rate.ddg.foldx')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voro   = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_elj, select = c(pdb,seqent,ddgent))
                  , subset(pdb_jec, select = c(zr4s_JC))
                  , subset(pdb_jec_ddg, select = c(rate.ddg.foldx))
                  , subset(pdb_hps, select = c(hpshh))
                  , subset(pdb_dssp, select = c(rsa,hbe))
                  , subset(pdb_wcn_bf, select = c(wcnSC,wcnCA,bfSC))
                  , subset(pdb_voro, select = c(VSCarea))
                  )
  r.rsa.seqent = cor(pdb_temp$seqent,pdb_temp$rsa,method='sp')
  r.wcnSC.seqent = cor(pdb_temp$seqent,pdb_temp$wcnSC,method='sp')
  r.wcnCA.seqent = cor(pdb_temp$seqent,pdb_temp$wcnCA,method='sp')
  r.vareaSC.seqent = cor(pdb_temp$seqent,pdb_temp$VSCarea,method='sp')
  #r.ddgent.seqent = cor(pdb_temp$seqent,pdb_temp$ddgent,method='sp')
  r.ddgent.seqent = cor(pdb_temp$seqent,pdb_temp$rate.ddg.foldx,method='sp')
  r.bfSC.seqent = cor(pdb_temp$seqent,pdb_temp$bfSC,method='sp')
  r.hbe.seqent = cor(pdb_temp$seqent,pdb_temp$hbe,method='sp')
  
  #row = data.frame( pdb = pdb, rsa = r.rsa.seqent, wcn = r.wcn.seqent, varea = r.vareaSC.seqent, ddgent = r.ddgent.seqent )
  row = data.frame( pdb, r.wcnSC.seqent, r.wcnCA.seqent, r.rsa.seqent, r.vareaSC.seqent, r.ddgent.seqent, r.bfSC.seqent )
  SEscors = rbind( SEscors, row )
}


# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
hist.vareaSC = density(SEscors$r.vareaSC.seqent)
hist.wcnSC = density(SEscors$r.wcnSC.seqent)
hist.wcnCA = density(SEscors$r.wcnCA.seqent)
hist.bfSC   = density(SEscors$r.bfSC.seqent)
hist.ddgent = density(SEscors$r.ddgent.seqent)
hist.rsa = density(SEscors$r.rsa.seqent)

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_structural_predictors_of_SE.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.vareaSC$x
       ,  hist.vareaSC$y
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
       ,   xlab = 'Spearman Cor. with Evolutionary Rates'
       ,   ylab = 'frequency'
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
        , c("Voronoi Cell Area (SC)","1 / WCN (SC)", "1 / WCN (CA)", "B factor (SC)", "ddG Rate", "RSA")
        , col = c('red','black','black','cyan2','green','blue')
        , lty = c(1,1,2,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)

graphics.off()
