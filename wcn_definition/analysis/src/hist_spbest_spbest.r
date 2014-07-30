# This R script uses data generated from my Fortran codes that I wrote in seach of the best performing definitions of WCN for individual structures compared to the original definitions, then creates histograms of the differences between the best performing correlations of the default-power-law definition of WCN and the alternative definitions.
# Inside the script, spbb stands for Best SPearman correlation vs. the spearman correlation of the default definition.

# Amir Shahmoradi, Friday 11:59 AM, July 18, 2014, Wilke Lab, iCMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_definition/analysis/src')   # Where this R code exists
# source('input_data.r')


####  GAUSSIAN WCN

    hist_spbb_wcng_bfac   = density(abs(sum_wcng_bfac$sp_best)-abs(sum_wcnp_bfac$sp_best))
    hist_spbb_wcng_seqent = density(abs(sum_wcng_seqent$sp_best)-abs(sum_wcnp_seqent$sp_best))
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_spbb_wcng.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_spbb_wcng_bfac$x, hist_spbb_wcng_bfac$y,
          col=colors[1],
          xlim = c(-.05,.05),
          #ylim = c(0.,50),
          #border = colors[1],
          #lty = 0,
          type = 'l',
          main = 'Gaussian Vs. Power-law WCN',
          xlab = 'diff. in Spearman cors',
          ylab = 'frequency'
    )
    # second histogram
    lines( hist_spbb_wcng_seqent$x, hist_spbb_wcng_seqent$y,
          col=rgb(1,0,0,1/4),
          xlim = c(-.05,.05),
          #ylim = c(0.,50),
          #border = colors[2],
          #lty = 0,
          #add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()
    

####  EXPONENTIAL WCN

    hist_spbb_wcne_bfac   = density(abs(sum_wcne_bfac$sp_best)-abs(sum_wcnp_bfac$sp_best))
    hist_spbb_wcne_seqent = density(abs(sum_wcne_seqent$sp_best)-abs(sum_wcnp_seqent$sp_best))
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_spbb_wcne.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_spbb_wcne_bfac$x, hist_spbb_wcne_bfac$y,
          col=colors[1],
          xlim = c(-.05,.05),
          #ylim = c(0.,50),
          #border = colors[1],
          #lty = 0,
          type = 'l',
          main = 'Exponential Vs. Power-law WCN',
          xlab = 'diff. in Spearman cors',
          ylab = 'frequency'
    )
    # second histogram
    lines( hist_spbb_wcne_seqent$x, hist_spbb_wcne_seqent$y,
          col=rgb(1,0,0,1/4),
          xlim = c(-.05,.05),
          #ylim = c(0.,50),
          #border = colors[2],
          #lty = 0,
          #add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()
    

####  DENSITY WCN

    hist_spbb_wcnd_bfac   = density(abs(sum_wcnd_bfac$sp_best)-abs(sum_wcnp_bfac$sp_best))
    hist_spbb_wcnd_seqent = density(abs(sum_wcnd_seqent$sp_best)-abs(sum_wcnp_seqent$sp_best))
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_spbb_wcnd.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_spbb_wcnd_bfac$x, hist_spbb_wcnd_bfac$y,
          col=colors[1],
          xlim = c(-.08,.08),
          #ylim = c(0.,50),
          #border = colors[1],
          #lty = 0,
          type = 'l',
          main = 'Density Vs. Power-law WCN',
          xlab = 'diff. in Spearman cors',
          ylab = 'frequency'
    )
    # second histogram
    lines( hist_spbb_wcnd_seqent$x, hist_spbb_wcnd_seqent$y,
          col=rgb(1,0,0,1/4),
          xlim = c(-.08,.08),
          #ylim = c(0.,50),
          #border = colors[2],
          #lty = 0,
          #add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()