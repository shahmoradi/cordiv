    
# NOW GENERATE CORRELATIONS HISTOGRAM DATA:
    hist_spcor.seqent_ddgent = density(pdb_prop_spcor$r.seqent_ddgent)
    hist_spcor.seqent_asa    = density(pdb_prop_spcor$r.seqent_asa)
    hist_spcor.seqent_rsa    = density(pdb_prop_spcor$r.seqent_rsa)
    hist_spcor.seqent_wcnca  = density(pdb_prop_spcor$r.seqent_wcnca)
    hist_spcor.seqent_hbe    = density(pdb_prop_spcor$r.seqent_hbe)
    hist_spcor.ddgent_asa    = density(pdb_prop_spcor$r.ddgent_asa)
    hist_spcor.ddgent_rsa    = density(pdb_prop_spcor$r.ddgent_rsa)
    hist_spcor.ddgent_wcnca  = density(pdb_prop_spcor$r.ddgent_wcnca)
    hist_spcor.ddgent_hbe    = density(pdb_prop_spcor$r.ddgent_hbe)
    hist_spcor.asa_rsa       = density(pdb_prop_spcor$r.asa_rsa)
    hist_spcor.asa_wcnca     = density(pdb_prop_spcor$r.asa_wcnca)
    hist_spcor.asa_hbe       = density(pdb_prop_spcor$r.asa_hbe)
    hist_spcor.rsa_wcnca     = density(pdb_prop_spcor$r.rsa_wcnca)
    hist_spcor.rsa_hbe       = density(pdb_prop_spcor$r.rsa_hbe)
    hist_spcor.wcnca_hbe     = density(pdb_prop_spcor$r.wcnca_hbe)
    
    #colors = c(rgb(0,0,1,1/4),rgb(1,0,0,1/4))
    colors = c('red', 'blue', 'green', 'purple', 'orange3', 'darkgreen', 'black', 'gray', 'cyan2')
    # Now plot both histograms in a single plot: first histogram
    pdf( "../figures/hist_spcor_seqent.pdf", width=4.5, height=4, useDingbats=FALSE )
    #par( mai=c(0.65, 0.65, 0.0, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
    
    # Add first histogram
    plot(  hist_spcor.seqent_ddgent$x,
           hist_spcor.seqent_ddgent$y,
           col=colors[1],
           xlim=c(-0.2,1),
           ylim=c(0,7),
           #border = colors[1],
           #lty = 0,
           type = 'l',
           lwd  = 2, 
           main = 'Correlations with Seqent',
           xlab = expression(paste('Spearman ',rho)),
           ylab = 'frequency'
    )
    # Add second histogram
    lines( hist_spcor.seqent_asa$x,
           hist_spcor.seqent_asa$y,
           col=colors[2],
           lwd  = 2, 
           #xlim = c(-.05,.05),
           #ylim = c(0.,50),
           #border = colors[2],
           #lty = 0,
           #add=T
           )
    # Add third histogram
    lines( hist_spcor.seqent_rsa$x,
           hist_spcor.seqent_rsa$y,
           col=colors[3],
           lwd  = 2, 
           #xlim = c(-.05,.05),
           #ylim = c(0.,50),
           #border = colors[2],
           #lty = 0,
           #add=T
           )
    # Add fourth histogram
    lines( hist_spcor.seqent_wcnca$x,
           hist_spcor.seqent_wcnca$y,
           col=colors[4],
           lwd  = 2, 
           #xlim = c(-.05,.05),
           #ylim = c(0.,50),
           #border = colors[2],
           #lty = 0,
           #add=T
           )
    # Add fifth histogram
    lines( hist_spcor.seqent_hbe$x,
           hist_spcor.seqent_hbe$y,
           col=colors[5],
           lwd  = 2, 
           #xlim = c(-.05,.05),
           #ylim = c(0.,50),
           #border = colors[2],
           #lty = 0,
           #add=T
           )
    legend('topright', c("WCN-Bfactor", "WCN-SeqEnt"), pch=19, col=colors[1:2], bty='n', cex=0.9 )
    #dev.off('all')
    #graphics.off()



dev.off()