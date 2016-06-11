# This R script creates the pipeline figure requested by the referee of the paper.
# Amir Shahmoradi, Tuesday 6:21 PM, January 26 2015, ICES, UT Austin

# Change directory to src directory:
setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# Source the script that reads the input data
source("get_res_data.r")

# 1. Generate a plot of ER-WCN and ER-Cell_volume for one single protein (here I am choosing the PDB with the strongest correlations, for the sake better illustration in paper)
pdb = "1ONR_A"
pdb_1ONR_A = data.frame( wcnSC = res_prop_wcn_bf$wcnSC[res_prop_wcn_bf$pdb == pdb]
                       , volumeSC = res_prop_voroSC$VSCvolume[res_prop_voroSC$pdb == pdb]
                       , r4sJC = res_prop_jec$zr4s_JC[res_prop_jec$pdb == pdb]
                       )

plot(log10(pdb_1ONR_A$wcnSC),log10(pdb_1ONR_A$volumeSC))

# Generate wcnSC-ER scatter plot for the PDB of interest
pdf( "../figures/pipeline/wcnSC-ER.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  pdb_1ONR_A$wcnSC
    ,  pdb_1ONR_A$r4sJC
    ,   pch  = 16
    ,   xlab = 'Weighted Contact Number (WCN)'
    ,   ylab = 'Site-Specific Evolutionary Rates (ER)'
)
graphics.off()

# Generate Cell_Volume-ER scatter plot for the PDB of interest
pdf( "../figures/pipeline/vvol-ER.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot(  log10(pdb_1ONR_A$volumeSC)
       ,  pdb_1ONR_A$r4sJC
       ,   pch  = 16
       ,   xlab = 'Log (Voronoi Cell Volume)'
       ,   ylab = 'Site-Specific Evolutionary Rates (ER)'
)
graphics.off()


# Now generate Cell_Volume-ER scatter plot for the PDB of interest
x = cor.test( pdb_1ONR_A$r4sJC, pdb_1ONR_A$wcnSC, method='spearman', na.action="na.omit" )
r.onr.wcn = x$estimate
x = cor.test( pdb_1ONR_A$r4sJC, pdb_1ONR_A$volumeSC, method='spearman', na.action="na.omit" )
r.onr.vvol = x$estimate
all_pdb_prop_select_wide = read.csv(file="../tables/all_pdb_prop_select_wide.csv", header=T)
pdf( "../figures/pipeline/ER_vvol_wcn.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot( all_pdb_prop_select_wide$r.r4sJC.vvolume
    , all_pdb_prop_select_wide$r.r4sJC.wcnSC
    , pch  = 16
    , xlab = expression( paste('Spearman ', rho , ':  ER-Cell Vol. Relation') )
    , ylab = expression( paste('Spearman ', rho , ':  ER-WCN Relation') )
    , xlim=c(0.25,0.85)
    , ylim=c(-0.25,-0.85)
)
points( r.onr.vvol
      , r.onr.wcn
      , col='red'
      , pch = 16
      )
abline(0,-1,col='red',lwd=2,lty=2)
graphics.off()


# Now generate Cell_Volume-ER scatter plot for the PDB of interest
#best_structural_predictors_of_ER_given_VSCvolume = read.csv(file = "../tables/best_structural_predictors_of_ER_given_VSCvolume.csv", header=T )
hist.ER.vvol = density(all_pdb_prop_select_wide$r.r4sJC.vvolume)
hist.ER.wcnSC = density(all_pdb_prop_select_wide$r.r4sJC.wcnSC)
pdf( "../figures/pipeline/ER-WCN-vvol_hist.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.65, 0.65, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot( -5
    , -5
    , col = 'red'
    , xlim = c(0.2,0.9)
    , ylim = c(0,5.)
    , type = 'l'
    , lwd  = 2 
    , xlab = expression( paste('Absolute Spearman ', rho ) )
    , ylab = 'Relative Frequency'
)
lines( -hist.ER.wcnSC$x
     , hist.ER.wcnSC$y
     , col = 'red'
     , lwd = 2
)
lines( hist.ER.vvol$x
     , hist.ER.vvol$y
     , col = 'black'
     , lwd = 2
)
legend( 'topleft'
        , c("ER - Cell Volume","ER - WCN")
        , col = c('black','red')
        , lty = c(1,1,2,1)
        , lwd = 2
        , bty = 'n'
        , cex = 0.9
)
graphics.off()



