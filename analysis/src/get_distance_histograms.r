# This R script generates histograms of the distributions of the nearest neighbor distances of sites in proteins.
# Amir Shahmoradi, Friday 9:32, May 1 2015, ICMB, UT Austin

# NOW GENERATE CORRELATIONS HISTOGRAM DATA:

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# First side chain distances

hist.nndSC = density(res_prop_nnd_SC$nnd)
hist.sndSC = density(res_prop_snd_SC$snd,na.rm=TRUE)
hist.nsnndSC = density(res_prop_nsnnd_SC$nsnnd)

pdf( "../figures/nearest_neighbor_distance_SC.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.sndSC$x
       #,  hist.sndSC$y
       ,  hist.sndSC$y*0.9/max(hist.sndSC$y)
       ,   col = 'black'
       #,   xlim = c(-0.5,0.85)
       ,   xlim = c(2.,12.5)
       ,   ylim=c(0,1.0)
       #,   ylim=c(0,0.015)
       #,   col=colors[1]
       #,   border = colors[1]
       #,   lty = 0
       ,   type = 'l'
       ,   lwd  = 2 
       #,   main = 'Correlations with Evolutionary Rates'
       #,   xlab = expression(paste('Spearman Cor. with Evolutionary Rates ',rho))
       ,   xlab = 'Distance [ Å ]'
       ,   ylab = 'Normalized Frequency'
)
lines(   hist.nsnndSC$x
       #, hist.nsnndSC$y
       , hist.nsnndSC$y*0.9/max(hist.nsnndSC$y)
       , col = 'blue'
       , lwd = 2
)
lines(   hist.nndSC$x
         #, hist.nndSC$y
         , hist.nndSC$y*0.9/max(hist.nndSC$y)
         , col = 'red'
         , lwd = 2
)
legend( 'topright'
        , c("Sequential Distance","NS-NN Distance","NN Distance")
        , col = c('black','blue','red')
        , lty = c(1,1,2,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)

graphics.off()


##############################
##############################
# Now do the same for CA atoms

hist.nndCA = density(res_prop_nnd_CA$nnd)
hist.sndCA = density(res_prop_snd_CA$snd,na.rm=TRUE)
hist.nsnndCA = density(res_prop_nsnnd_CA$nsnnd)

pdf( "../figures/nearest_neighbor_distance_CA.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  hist.sndCA$x
       ,  hist.sndCA$y*0.9/max(hist.sndCA$y)
       ,   col = 'black'
       #,   xlim = c(3.5,7.5)
       ,   xlim = c(2.,12.5)
       ,   ylim=c(0.0,1.0)
       #,   col=colors[1]
       #,   border = colors[1]
       #,   lty = 0
       ,   type = 'l'
       ,   lwd  = 2 
       #,   main = 'Correlations with Evolutionary Rates'
       #,   xlab = expression(paste('Spearman Cor. with Evolutionary Rates ',rho))
       ,   xlab = 'Distance [ Å ]'
       ,   ylab = 'Normalized Frequency'
       )
lines(   hist.nsnndCA$x
         , hist.nsnndCA$y*0.9/max(hist.nsnndCA$y)
         , col = 'blue'
         , lwd = 2
)
lines(   hist.nndCA$x
         , hist.nndCA$y*0.9/max(hist.nndCA$y)
         , col = 'red'
         , lwd = 2
)
legend( 'topright'
        , c("Sequential Distance","NS-NN Distance","NN Distance")
        , col = c('black','blue','red')
        , lty = c(1,1,2,1,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)

graphics.off()
