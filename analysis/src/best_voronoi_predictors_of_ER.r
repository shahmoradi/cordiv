# This code aims at finding which Voronoi cell property best predicts the Evolutionary Rates.

# Last updated by Amir Shahmoradi, Thursday 9:38 AM, February 26 2015, Wilke Lab, ICMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
library(ppcor)
pdb_temp = data.frame()
best_voronoi_predictors_of_ER = data.frame()  # This data frame will contain correlations of select structural variables with r4sJC evolutionary rates, for each pdb file on eaxch row
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  #pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
  #pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  #pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
  #pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  #pdb_voroCA = res_prop_voroCA[res_prop_voroCA$pdb==pdb, ]
  
  pdb_temp = cbind( subset(pdb_jec, select = c(zr4s_JC))
                  #, subset(pdb_elj, select = c(pdb,seqent,ddgent))
                  #, subset(pdb_hps, select = c(hpshh))
                  #, subset(pdb_dssp, select = c(rsa,hbe))
                  #, subset(pdb_wcn_bf, select = c(wcnSC,bfSC))
                  , subset(pdb_voroSC, select = c(VSCsphericity,VSCeccentricity,VSCvolume,VSCarea,VSCedge_length_total,VSCnfaces,VSCnedges,VSCnvertices))
                  )
  #nvertices = cor(pdb_temp$zr4s_JC,pdb_temp$VSCnvertices,method='sp') # This is identical to the two folloing variables and is therfore dropped from the list
  #nedges = cor(pdb_temp$zr4s_JC,pdb_temp$VSCnedges,method='sp')
  #nfaces = cor(pdb_temp$zr4s_JC,pdb_temp$VSCnfaces,method='sp')
  edge = cor(pdb_temp$zr4s_JC,pdb_temp$VSCedge_length_total,method='sp')
  area = cor(pdb_temp$zr4s_JC,pdb_temp$VSCarea,method='sp')
  volume = cor(pdb_temp$zr4s_JC,pdb_temp$VSCvolume,method='sp')
  eccentricity = cor(pdb_temp$zr4s_JC,pdb_temp$VSCeccentricity,method='sp')
  sphericity = cor(pdb_temp$zr4s_JC,pdb_temp$VSCsphericity,method='sp')
  
  # out of curiosity: test for partial correlation of cell variables with ER, controlling for cell area:
  x = pcor.test( pdb_temp$zr4s_JC, pdb_temp$VSCedge, pdb_temp$VSCarea, method='sp')
  r.edg_given_area = x$estimate
  x = pcor.test( pdb_temp$zr4s_JC, pdb_temp$VSCvolume, pdb_temp$VSCarea, method='sp')
  r.vol_given_area = x$estimate
  x = pcor.test( pdb_temp$zr4s_JC, pdb_temp$VSCeccentricity, pdb_temp$VSCarea, method='sp')
  r.ecc_given_area = x$estimate
  x = pcor.test( pdb_temp$zr4s_JC, pdb_temp$VSCsphericity, pdb_temp$VSCarea, method='sp')
  r.sph_given_area = x$estimate
  #pcor.test( pdb_temp$zr4s_JC, pdb_temp$VSCeccentricity, pdb_temp$VSCsphericity, method='sp')
  #row = data.frame( pdb = pdb, rsa = r.rsa.r4sJC, wcn = r.wcn.r4sJC, vareaSC = r.vareaSC.r4sJC, ddgent = r.ddgent.r4sJC )
  #row = data.frame( pdb, r.nvertices.r4sJC, r.nedges.r4sJC, r.nfaces.r4sJC, r.edgelength.r4sJC, r.area.r4sJC , r.volume.r4sJC, r.eccentricity.r4sJC, r.sphericity.r4sJC )
  row = data.frame( pdb, edge, area , volume, eccentricity, sphericity
                  , r.edg_given_area, r.vol_given_area, r.ecc_given_area, r.sph_given_area )
  best_voronoi_predictors_of_ER = rbind( best_voronoi_predictors_of_ER, row )
}

write.csv(best_voronoi_predictors_of_ER, file = "../tables/best_voronoi_predictors_of_ER.csv", row.names=F )

#GGPLOT is a waste of time!
### library(reshape2)
### library(ggplot2)
# The palette with black:
cbbPalette <- c("#009E73", "#0072B2", "#D55E00", "#CC79A7", "#000000", "#F0E442", "#E69F00") #, "#56B4E9")
# To use for fills, add
#scale_fill_manual(values=cbPalette)
# To use for line and point colors, add
#scale_colour_manual(values=cbPalette)

### best_voronoi_predictors_of_ER_long = melt(best_voronoi_predictors_of_ER, id.vars = "pdb")
### pdf( "../figures/best_voronoi_predictors_of_ER.pdf", width=4.5, height=4, useDingbats=FALSE )
### par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
### ggplot(best_voronoi_predictors_of_ER_long, aes(x=value, colour=variable) ) + geom_density(size = 1, fill=NA) + scale_colour_manual(values=cbbPalette) + xlim(0.1,0.9) + 
### graphics.off()


# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
hist.edge = density(best_voronoi_predictors_of_ER$edge)
hist.area = density(best_voronoi_predictors_of_ER$area)
hist.volume = density(best_voronoi_predictors_of_ER$volume)
hist.eccentricity = density(best_voronoi_predictors_of_ER$eccentricity)
hist.sphericity = density(best_voronoi_predictors_of_ER$sphericity)
hist.r.edg_given_area = density(best_voronoi_predictors_of_ER$r.edg_given_area)
hist.r.vol_given_area = density(best_voronoi_predictors_of_ER$r.vol_given_area)
hist.r.ecc_given_area = density(best_voronoi_predictors_of_ER$r.ecc_given_area)
hist.r.sph_given_area = density(best_voronoi_predictors_of_ER$r.sph_given_area)

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_voronoi_predictors_of_ER.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.edge$x
    ,  hist.edge$y
    #, col = 'blue'
    ,  col = cbbPalette[[1]]
    ,   xlim = c(0.1,0.85)
    ,   ylim = c(0,5.5)
    #,   col=colors[1]
    #,   ylim=c(0,7)
    #,   border = colors[1]
    #,   lty = 0
    ,   type = 'l'
    ,   lwd  = 2 
    #,   main = 'Correlations with Evolutionary Rates'
    #,   xlab = expression(paste('Absolute Spearman Cor. with Evolutionary Rates ',rho))
    ,   xlab = 'Absolute Spearman Cor. with Evolutionary Rates'
    ,   ylab = 'frequency'
    )
lines( abs(hist.area$x)
     , abs(hist.area$y)
     #, col = 'green'
     ,  col = cbbPalette[[2]]
     , lwd  = 2
     )
lines( abs(hist.volume$x)
     , abs(hist.volume$y)
     #, col = 'red'
     , col = cbbPalette[[3]]
     , lwd  = 2
     )
lines( abs(hist.eccentricity$x)
     , abs(hist.eccentricity$y)
     #, col = 'black'
     , col = cbbPalette[[4]]
     , lwd = 2
     )
lines( abs(hist.sphericity$x)
     , abs(hist.sphericity$y)
     #, col = 'black'
     , col = cbbPalette[[5]]
     , lwd = 2
     #, lty = 2
     )
legend( 'topleft'
        , c("cell edge length", "cell area", "cell volume", "cell eccentricity", "cell sphericity")
        #, col = c('red','black','black','green','blue')
        , col = cbbPalette
        , lty = c(1,1,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
      )


# Now plot histograms of correlations with ER, while controlling for cell area:
# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_voronoi_predictors_of_ER_given_area.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.r.edg_given_area$x
    ,  hist.r.edg_given_area$y
    #, col = 'blue'
    ,  col = cbbPalette[[1]]
    ,   xlim = c(-0.4,0.3)
    ,   ylim = c(0,7.5)
    #,   col=colors[1]
    #,   ylim=c(0,7)
    #,   border = colors[1]
    #,   lty = 0
    ,   type = 'l'
    ,   lwd  = 2 
    #,   main = 'Correlations with Evolutionary Rates'
    #,   xlab = expression(paste('Absolute Spearman Cor. with Evolutionary Rates ',rho))
    ,   xlab = 'Spearman Cor. with Evolutionary Rates'
    ,   ylab = 'frequency'
    )
lines( hist.r.vol_given_area$x
     , hist.r.vol_given_area$y
     #, col = 'red'
     , col = cbbPalette[[2]]
     , lwd  = 2
     )
lines( hist.r.ecc_given_area$x
     , hist.r.ecc_given_area$y
     #, col = 'black'
     , col = cbbPalette[[3]]
     , lwd = 2
     #, lty = 2
     )
lines( hist.r.sph_given_area$x
       , hist.r.sph_given_area$y
       #, col = 'black'
       , col = cbbPalette[[4]]
       , lwd = 2
       #, lty = 2
     )
legend( 'topleft'
      , c("cell edge length", "cell volume", "cell eccentricity", "cell sphericity")
      #, col = c('red','black','green','blue')
      , col = cbbPalette
      , lty = c(1,1,1)
      , lwd = 2
      , bty = 'n'
      , cex = 0.9
      )

graphics.off()




#####################################################
#####################################################
#####################################################
#####################################################
#####################################################



# Now combine both histograms in one screen figure:

# Now plot histograms in a single plot
#colors = c('green', 'blue', 'red', 'black', 'gray', 'cyan2')
pdf( "../figures/best_voronoi_predictors_of_ER_screen.pdf", width=9, height=4, useDingbats=FALSE )

split.screen(c(1,2))

  screen(1)
  par( mai=c(0.65, 0.65, 0.2, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot( hist.edge$x
    , hist.edge$y
    #, col = 'blue'
    ,  col = cbbPalette[[1]]
    ,   xlim = c(0.1,0.85)
    ,   ylim = c(0,5.5)
    #,   col=colors[1]
    #,   ylim=c(0,7)
    #,   border = colors[1]
    #,   lty = 0
    ,   type = 'l'
    ,   lwd  = 2 
    #,   main = 'Correlations with Evolutionary Rates'
    #,   xlab = expression(paste('Absolute Spearman Cor. with Evolutionary Rates ',rho))
    ,   xlab = 'Absolute Spearman Cor. with Evolutionary Rates'
    ,   ylab = 'frequency'
    )
  mtext('A', side = 3, at=-0.03, font=2, cex=1.2)
  lines( abs(hist.area$x)
       , abs(hist.area$y)
       #, col = 'green'
       ,  col = cbbPalette[[2]]
       , lwd  = 2
      )
  lines( abs(hist.volume$x)
       , abs(hist.volume$y)
       #, col = 'red'
       , col = cbbPalette[[3]]
       , lwd  = 2
       )
  lines( abs(hist.eccentricity$x)
       , abs(hist.eccentricity$y)
       #, col = 'black'
       , col = cbbPalette[[4]]
       , lwd = 2
       )
  lines( abs(hist.sphericity$x)
       , abs(hist.sphericity$y)
       #, col = 'black'
       , col = cbbPalette[[5]]
       , lwd = 2
       #, lty = 2
       )
  legend( 'topleft'
        , c("cell edge length", "cell area", "cell volume", "cell eccentricity", "cell sphericity")
        #, col = c('red','black','black','green','blue')
        , col = cbbPalette
        , lty = c(1,1,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
        )


screen(2)

  par( mai=c(0.65, 0.65, 0.2, 0.2), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(  hist.r.edg_given_area$x
       ,  hist.r.edg_given_area$y
       #, col = 'blue'
       ,  col = cbbPalette[[1]]
       ,   xlim = c(-0.4,0.3)
       ,   ylim = c(0,7.5)
       #,   col=colors[1]
       #,   ylim=c(0,7)
       #,   border = colors[1]
       #,   lty = 0
       ,   type = 'l'
       ,   lwd  = 2 
       #,   main = 'Correlations with Evolutionary Rates'
       #,   xlab = expression(paste('Absolute Spearman Cor. with Evolutionary Rates ',rho))
       ,   xlab = 'Spearman Cor. with Evolutionary Rates'
       ,   ylab = 'frequency'
       )
  mtext('B', side = 3, at=-0.52, font=2, cex=1.2)
  lines( hist.r.vol_given_area$x
       , hist.r.vol_given_area$y
       #, col = 'red'
       , col = cbbPalette[[2]]
       , lwd  = 2
       )
  lines( hist.r.ecc_given_area$x
       , hist.r.ecc_given_area$y
       #, col = 'black'
       , col = cbbPalette[[3]]
       , lwd = 2
       #, lty = 2
       )
  lines( hist.r.sph_given_area$x
       , hist.r.sph_given_area$y
       #, col = 'black'
       , col = cbbPalette[[4]]
       , lwd = 2
       #, lty = 2
       )
  legend( 'topleft'
        , c("cell edge length", "cell volume", "cell eccentricity", "cell sphericity")
        #, col = c('red','black','green','blue')
        , col = cbbPalette
        , lty = c(1,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
        )

close.screen(all = TRUE)

graphics.off()

