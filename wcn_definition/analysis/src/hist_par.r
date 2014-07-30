# This R script uses data generated from my Fortran codes that I wrote in seach of the best performing definitions of WCN for individual structures compared to the original definitions, then creates histograms of the best performing parameters and correlations.
# Amir Shahmoradi, Friday 10:35 AM, July 18, 2014, Wilke Lab, iCMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_definition/analysis/src')   # Where this R code exists
# source('input_data.r')


#### POWER LAW WCN

    hist_par_wcnp_bfac   = hist(sum_wcnp_bfac$exp_best, breaks=200)
    hist_par_wcnp_seqent = hist(sum_wcnp_seqent$exp_best, breaks=200)
    
    # Now plot both histograms in a single plot: first histogram
      pdf( "../figures/hist_par_wcnp.pdf", width=4.5, height=4, useDingbats=FALSE )
      #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
      colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
      plot( hist_par_wcnp_bfac,
            col=colors[1],
            xlim = c(-5.,30),
            ylim = c(0.,25),
            #border = colors[1],
            lty = 0,
            main = 'Power-law WCN',
            xlab = 'best performing exponent'
            )
    # second histogram
      plot( hist_par_wcnp_seqent,
            col=rgb(1,0,0,1/4),
            xlim = c(-5.,30),
            ylim = c(0.,25),
            #border = colors[2],
            lty = 0,
            add=T
            )
      legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
      #dev.off('all')
      graphics.off()
    
    # Now plot both histograms zoomed in a narrow range for better visibility, in a single plot again:
    pdf( "../figures/hist_par_wcnp_zoomedin.pdf", width=4.5, height=4, useDingbats=FALSE )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcnp_bfac,
          col=colors[1],
          xlim = c(-5.,2.),
          ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'Power-law WCN',
          xlab = 'best performing exponent'
    )
    # second histogram
    plot( hist_par_wcnp_seqent,
          col=rgb(1,0,0,1/4),
          xlim = c(-5.,2.),
          ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()


####  GAUSSIAN WCN

    hist_par_wcng_bfac   = hist(sum_wcng_bfac$std_best, breaks=50)
    hist_par_wcng_seqent = hist(sum_wcng_seqent$std_best, breaks=50)
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_par_wcng.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcng_bfac,
          col=colors[1],
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'Gaussian WCN',
          xlab = 'best performing distance [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcng_seqent,
          col=rgb(1,0,0,1/4),
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()
    
    # Now plot both histograms zoomed in a narrow range for better visibility, in a single plot again:
    pdf( "../figures/hist_par_wcng_zoomedin.pdf", width=4.5, height=4, useDingbats=FALSE )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcng_bfac,
          col=colors[1],
          xlim = c(0.,25.),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'Gaussian WCN',
          xlab = 'best performing distance [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcng_seqent,
          col=rgb(1,0,0,1/4),
          xlim = c(0.,25.),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()


####  EXPONENTIAL WCN

    hist_par_wcne_bfac   = hist(sum_wcne_bfac$expmean_best, breaks=50)
    hist_par_wcne_seqent = hist(sum_wcne_seqent$expmean_best, breaks=50)
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_par_wcne.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcne_bfac,
          col=colors[1],
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'Exponential WCN',
          xlab = 'best performing mean [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcne_seqent,
          col=rgb(1,0,0,1/4),
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()
    
    # Now plot both histograms zoomed in a narrow range for better visibility, in a single plot again:
    pdf( "../figures/hist_par_wcne_zoomedin.pdf", width=4.5, height=4, useDingbats=FALSE )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcne_bfac,
          col=colors[1],
          xlim = c(0.,15.),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'Exponential WCN',
          xlab = 'best performing mean [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcne_seqent,
          col=rgb(1,0,0,1/4),
          xlim = c(0.,15.),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()


####  DENSITY WCN

    hist_par_wcnd_bfac   = hist(sum_wcnd_bfac$cutoff_best, breaks=50)
    hist_par_wcnd_seqent = hist(sum_wcnd_seqent$cutoff_best, breaks=50)
    
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_par_wcnd.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcnd_bfac,
          col=colors[1],
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'cutoff CN',
          xlab = 'best performing cutoff [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcnd_seqent,
          col=rgb(1,0,0,1/4),
          #xlim = c(-5.,30),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()
    
    # Now plot both histograms zoomed in a narrow range for better visibility, in a single plot again:
    pdf( "../figures/hist_par_wcnd_zoomedin.pdf", width=4.5, height=4, useDingbats=FALSE )
    colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    plot( hist_par_wcnd_bfac,
          col=colors[1],
          xlim = c(5.,30.),
          #ylim = c(0.,25),
          #border = colors[1],
          lty = 0,
          main = 'cutoff CN',
          xlab = 'best performing cutoff [Angstroms]'
    )
    # second histogram
    plot( hist_par_wcnd_seqent,
          col=rgb(1,0,0,1/4),
          xlim = c(5.,30.),
          #ylim = c(0.,25),
          #border = colors[2],
          lty = 0,
          add=T
    )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    graphics.off()














# wcn_bfac_gaus = read.table('../wcn/gaussian/bfactor/sum_wcng_bfac.out', header=T)
# hist(wcn_bfac_gaus$std_best, breaks = 25, xlim = c(0,50))
# hist(wcn_bfac_gaus$sp_diff,breaks=20)
# plot(abs(wcn_bfac_gaus$sp_best),abs(wcn_bfac_pwrl$sp_best))
# 
# wcn_seqent_gaus = read.table('../wcn/gaussian/seqent/sum_wcng_seqent.out', header=T)
# hist_std_seqent = hist(wcn_seqent_gaus$std_best, breaks = 25, xlim = c(0,50))
# hist(wcn_seqent_gaus$sp_diff,breaks=20)
# hist(abs(wcn_seqent_gaus$sp_best)-abs(wcn_seqent_pwrl$sp_best),breaks=20)
# 
# wcn_bfac_expn = read.table('../wcn/exponential/bfactor/sum_wcne_bfac.out', header=T)
# hist(wcn_bfac_expn$expmean_best, breaks = 25, xlim = c(0,50))
# hist(wcn_bfac_expn$sp_diff,breaks=20)
# mean(wcn_bfac_expn$sp_diff)
# hist(abs(wcn_bfac_expn$sp_best)-abs(wcn_bfac_pwrl$sp_best),breaks=20,xlim = c(-.15,0.05))
# median(abs(wcn_bfac_expn$sp_best)-abs(wcn_bfac_pwrl$sp_best))
# mean(abs(wcn_bfac_expn$sp_best)-abs(wcn_bfac_gaus$sp_best))