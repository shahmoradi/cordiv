# This R code generates the plots of correlations of variables with seqent, vs. sd.seqent and compare it to ASAP data.
# Amir Shahmoradi, Thursday 6:38 PM, Oct 9, 2014, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# source('input_ASAP_data')

library('ppcor')

ASAP_pdb_data = read.csv("../tables/ASAP_pdb_prop.csv", header = T)

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv", header = T)
#View(all_pdb_prop_select_wide)

n.screen.rows = 2   # number of rows in the screen plots
n.screen.cols = 2   # number of columns in the screen plots

# First generate a screen plot of the correlation of seqent-wcnSC vs.seqent-other

filename = paste0('../figures/cordiv_similarities_seqent.pdf')
pdf( filename, width=4.5*n.screen.cols, height=4*n.screen.rows, useDingbats=FALSE )
column_list = c('r.rsa.seqent'
                ,'r.ddgent.seqent'
                ,'r.bfSC.seqent'
                ,'r.seqent.varea'
                ,'r.seqent.vsphericity'
                ,'r.seqent.veccentricity'
                )
name_list = c('( Seq. Entropy - RSA )'
              ,'( Seq. Entropy - ddG Entropy )'
              ,'( Seq. Entropy - Mean Residue Bfactor )'
              ,'( Seq. Entropy - Voronoi Cell Surface Area )'
              ,'( Seq. Entropy - Voronoi Cell Compactness )'
              ,'( Seq. Entropy - Voronoi Cell Symmetry )'
              )
split.screen(c(n.screen.rows,n.screen.cols))
x = -1:1
for (counter in seq(1,n.screen.rows*n.screen.cols))
{
  screen(counter)
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(abs(all_pdb_prop_select_wide[[column_list[counter]]]),
       abs(all_pdb_prop_select_wide$r.seqent.wcnSC),
       xlab = bquote('abs.' ~ rho ~ .(name_list[counter])),
       ylab = expression(paste('abs. ', rho,' ( Seq. Entropy - Side Chain WCN )')),
       xlim = c(0.,0.9),
       ylim = c(0.,0.9),
       #cex.axis = 1.4,
       #cex.lab = 1.4,
       pch = 16
       )
  lines(x,x,col='red')
}
dev.off()
#graphics.off()

# Now generate a similar screen plot of the correlation of r4sJC-wcnSC vs. r4sJC-other
filename = paste0('../figures/cordiv_similarities_r4sJC.pdf')
pdf( filename, width=4.5*n.screen.cols, height=4*n.screen.rows, useDingbats=FALSE )
column_list = c('r.r4sJC.varea'
               ,'r.r4sJC.rsa'
               ,'r.ddgent.r4sJC'
               ,'r.bfSC.r4sJC'
               ,'r.r4sJC.vsphericity'
               ,'r.r4sJC.veccentricity'
               )
name_list = c('( Evol. Rates - Voronoi Cell Surface Area )'
              ,'( Evol. Rates - RSA )'
              ,'( Evol. Rates - ddG Entropy )'
              ,'( Evol. Rates - Mean Residue Bfactor )'
              ,'( Evol. Rates - Voronoi Cell Compactness )'
              ,'( Evol. Rates - Voronoi Cell Symmetry )'
              )
split.screen(c(n.screen.rows,n.screen.cols))
x = -1:1
for (counter in seq(1,n.screen.rows*n.screen.cols))
{
  screen(counter)
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(abs(all_pdb_prop_select_wide[[column_list[counter]]]),
       abs(all_pdb_prop_select_wide$r.r4sJC.wcnSC),
       xlab = bquote('abs.' ~ rho ~ .(name_list[counter])),
       ylab = expression(paste('abs. ', rho,' ( Evol. Rates - Side Chain WCN )')),
       xlim = c(0.,0.9),
       ylim = c(0.,0.9),
       #cex.axis = 1.4,
       #cex.lab = 1.4,
       pch = 16
       )
  lines(x,x,col='red')
  scor = cor(abs(all_pdb_prop_select_wide[[column_list[counter]]]),abs(all_pdb_prop_select_wide$r.r4sJC.wcnSC))
  text(x = 0.70, y = 0.05, labels = bquote('Spearman' ~ rho ~ '=' ~ .(round(scor,2))))
}
dev.off()
#graphics.off()




# Then generate a screen plot of the correlation of sd.seqent with seqent-variable correlations
counter = 0
filename = paste0('../figures/sd.seqent_cors.pdf')
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


# Now plot data for validation

source('input_ASAP_data.r')

filename = paste0('../figures/validation_sdseqent.pdf')
pdf( filename, width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )

plot( all_pdb_prop_select_wide$sd.seqent
    , all_pdb_prop_select_wide$r.seqent.wcnSC
    , pch = 19
    , xlab = 'St. Dev. of Seq. Entropy ( sd.seqent )'
    , ylab = bquote('Spearman' ~ rho ~ 'seqent - wcnSC')
    , xlim = c(-0.05,0.9)
    , ylim = c(-0.8,0.1)
    )
points( ASAP_pdb_prop$sd.seqent
      , ASAP_pdb_prop$r.seqent.wcnCA
      , col = 'red'
      , pch = 19
      )
dev.off()

filename = paste0('../figures/validation_sdhbe.pdf')
pdf( filename, width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )

plot( all_pdb_prop_select_wide$sd.hbe
      , all_pdb_prop_select_wide$r.r4sJC.wcnSC
      , pch = 19
      , xlab = 'St. Dev. of Seq. Entropy ( sd.seqent )'
      , ylab = bquote('Spearman' ~ rho ~ 'seqent - wcnSC')
      #, xlim = c(-0.05,0.9)
      #, ylim = c(-0.8,0.1)
)
# points( ASAP_pdb_prop$sd.seqent
#         , ASAP_pdb_prop$r.seqent.wcnCA
#         , col = 'red'
#         , pch = 19
# )
dev.off()


# For Claus proposal:

source('input_ASAP_data.r')
filename = paste0('../figures/cordiv_seqent_wcn_rsa.pdf')
pdf( filename, width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.1, 0.1), mgp=c(2, 0.5, 0), tck=-0.03 )
plot( all_pdb_prop_select_wide$r.rsa.seqent
    , all_pdb_prop_select_wide$r.seqent.wcnSC
    , pch = 19
    , xlab = bquote('Spearman' ~ rho ~ '( seq. entropy - RSA )')
    , ylab = bquote('Spearman' ~ rho ~ '( seq. entropy - WCN )')
    , xlim = c(-0.1,0.85)
    , ylim = c(0.1,-0.85)
)
x = -1:1
lines(x,-x,col='green',lwd=2)
points( ASAP_pdb_prop$r.seqent.rsa
      , ASAP_pdb_prop$r.seqent.wcnCA
      , pch = 19
      , col = 'red'
      )
dev.off()

plot (abs(all_pdb_prop_select_wide$r.seqent.wcnSC), abs(all_pdb_prop_select_wide$r.seqent.vsphericitym))
lines(x,x)
plot (abs(all_pdb_prop_select_wide$r.seqent.wcnSC), abs(all_pdb_prop_select_wide$r.seqent.vsphericitym))
lines(x,x)
plot (abs(all_pdb_prop_select_wide$r.seqent.wcnSC), abs(all_pdb_prop_select_wide$r.rsa.seqent))
lines(x,x)
plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.varea)
plot (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.seqent.veccentricity)
plot (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.ddgent.seqent)
plot (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent)
plot (all_pdb_prop_select_wide$sd.vccsphericity, all_pdb_prop_select_wide$mean.vccsphericity)
cor.test (all_pdb_prop_select_wide$mean.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, method='sp')
cor.test (all_pdb_prop_select_wide$sd.seqent, all_pdb_prop_select_wide$r.rsa.seqent, method='sp')
pcor.test (all_pdb_prop_select_wide$mean.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, all_pdb_prop_select_wide$sd.mvsphericity, method='sp')
pcor.test (all_pdb_prop_select_wide$sd.mvsphericity, all_pdb_prop_select_wide$r.mvsphericity.seqent, all_pdb_prop_select_wide$mean.mvsphericity, method='sp')
cor.test (all_pdb_prop_select_wide$contact_orderSC, all_pdb_prop_select_wide$r.ddgent.seqent)
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
