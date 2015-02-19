# This code is aimed at finding which structural property best predicts  and area (varea) can best represent a residue. This is done by comparing the correlations of different vorvols with other residue variables.

# Last updated by Amir Shahmoradi, Friday 4:25 PM, February 13 2015, Wilke Lab, ICMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
all_scors_select_variables = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  if (!(pdb %in% excluded_pdbs))
  {
    counter = counter + 1
    cat( paste(str(counter),pdb) )
    
    pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
    pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
    pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
    pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
    pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ] #c('asa','rsa','hbe_mean','rss')] )
    
    pdb_temp = cbind( subset(pdb_elj, select = c(seqent,ddgent)),
                      subset(pdb_jec, select = c(zr4s_JC)),
                      subset(pdb_hps, select = c(hpskd,hpsww,hpshh)),
                      subset(pdb_dssp, select = c(asa,rsa,hbe_mean)),
                      subset(pdb_wcn_bf, select = -c(pdb,resnam,resnum))
    )
    
    pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
    pdb_long$variable = factor(pdb_long$variable)
    
    counter1 = 0
    
    for (variable1 in levels(pdb_long$variable))
    {
      counter1 = counter1 + 1
      #cat (variable1, '\n')
      var1 = pdb_long[pdb_long$variable == variable1,]
      
      # Claculate potentially important statistical moments of the factored variable:
      row = data.frame(pdb, variable = paste0('sum.',variable1), value = sum(var1$value))       ; pdb_prop_from_residue_prop = rbind(pdb_prop_from_residue_prop,row)
      row = data.frame(pdb, variable = paste0('mean.',variable1), value = mean(var1$value))     ; pdb_prop_from_residue_prop = rbind(pdb_prop_from_residue_prop,row)
      row = data.frame(pdb, variable = paste0('median.',variable1), value = median(var1$value)) ; pdb_prop_from_residue_prop = rbind(pdb_prop_from_residue_prop,row)
      row = data.frame(pdb, variable = paste0('sd.',variable1), value = sd(var1$value))         ; pdb_prop_from_residue_prop = rbind(pdb_prop_from_residue_prop,row)
      
      # Now calculate the Spearman correlations between pairs of variables:
      counter2 = 0
      for (variable2 in levels(pdb_long$variable))
      {
        counter2 = counter2 + 1
        if ( variable1 != variable2 & counter1 < counter2)
        {
          var2 = pdb_long[pdb_long$variable == variable2,]
          x = cor.test( var1$value, var2$value, method='spearman', na.action="na.omit" )
          r = x$estimate
          p = x$p.value
          
          row = data.frame(pdb, variable = paste0('r.',variable1,'.',variable2), value = r)
          pdb_prop_from_residue_prop = rbind(pdb_prop_from_residue_prop,row)
        }
      }
    }
  }
}


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