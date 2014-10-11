# This R code generates the plots of correlations of variables with seqent, vs. sd.seqent and compare it to ASAP data.
# Amir Shahmoradi, Thursday 6:38 PM, Oct 9, 2014, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# source('input_ASAP_data')

library('ppcor')

ASAP_pdb_data = read.csv("../tables/ASAP_pdb_prop.csv", header = T)

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv", header = T)
#View(all_pdb_prop_select_wide)

# First generate a screen plot of the correlation of seqent-wcnSC vs.seqent-other
counter = 0
filename = paste0('../figures/seqent_structure_cors.pdf')
pdf( filename, width=9, height=8, useDingbats=FALSE )
column_list = c('r.rsa.seqent','r.ddgent.seqent','r.bfSC.seqent','r.seqent.varea')
name_list = c('rho ( Seq. Entropy - RSA )','rho ( Seq. Entropy - ddG Entropy )','rho ( Seq. Entropy - Bfactor )', 'rho ( Seq. Entropy - Voronoi Cell Volume )')
split.screen(c(2,2))
x = -1:1
for (column in column_list)
{
  counter = counter + 1
  screen(counter)
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(all_pdb_prop_select_wide[[column_list[counter]]],
       abs(all_pdb_prop_select_wide$r.seqent.wcnSC),
       xlab = name_list[counter],
       ylab = 'absolute rho ( Seq. Entropy - wcnSC )',
       xlim = c(-0.1,0.9),
       ylim = c(-0.1,0.9)
       )
  lines(x,x,col='red')
}
dev.off()
#graphics.off()

plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.varea)
plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.veccentricity)
plot (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.ddgent.seqent)
plot (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent)
plot (all_pdb_prop_select_wide$sd.vccsphericity, all_pdb_prop_select_wide$mean.vccsphericity)
cor.test (all_pdb_prop_select_wide$mean.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, method='sp')
cor.test (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.rsa.seqent, method='sp')
pcor.test (all_pdb_prop_select_wide$mean.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, all_pdb_prop_select_wide$sd.mvsphericity, method='sp')
pcor.test (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, all_pdb_prop_select_wide$mean.mvsphericity, method='sp')
cor.test (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.ddgent.seqent)
cor.test (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.ddgent.mvsphericity)
plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.mvsphericity.seqent)
cor.test (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.veccentricity, method='sp')
plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.ddgent.seqent)

plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.rsa.seqent, xlim = c(0,0.9), ylim = c(-0.1,1))
points(ASAP_pdb_data$sd.seqent, ASAP_pdb_data$r.seqent.rsa)

plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.bfSC.seqent, xlim = c(0,0.8), ylim = c(-0.2,1))
points(ASAP_pdb_data$sd.seqent, ASAP_pdb_data$r.seqent.bfCA)

plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.wcnSC, xlim = c(0,1.2), ylim = c(-1,0))
points(ASAP_pdb_data$sd.seqent, ASAP_pdb_data$r.seqent.wcnCA)

cor.test (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$mean.seqent)
plot (all_pdb_prop_select_wide$sum.seqent, all_pdb_prop_select_wide$r.seqent.varea)
plot (all_pdb_prop_select_wide$mean.seqent, all_pdb_prop_select_wide$r.seqent.varea)

pcor.test(all_pdb_prop_select_wide$mean.seqent, all_pdb_prop_select_wide$r.seqent.wcnSC, all_pdb_prop_select_wide$sd.seqent, method='sp')
pcor.test(all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.wcnSC, all_pdb_prop_select_wide$mean.seqent, method='sp')
cor.test(all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.wcnSC, method='sp')
cor.test(all_pdb_prop_select_wide$mean.seqent, all_pdb_prop_select_wide$r.seqent.wcnSC, method='sp')
cor.test(all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$sum.seqent, method='sp')
