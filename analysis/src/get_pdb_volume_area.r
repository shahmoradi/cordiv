# Last updated by Amir Shahmoradi, Friday 5:59 PM, Feb 20 2015, Wilke Lab, ICMB, UT Austin

# This R script reads in voluems and surface areas of pdbs in the dataset of 209 proteins. the volume and areas are calculated using 3V software written by,
#   Voss NR, Gerstein M, Steitz TA, Moore PB.
#   "The geometry of the ribosomal polypeptide exit tunnel."
#   J Mol Biol. 2006 Jul 21;360(4):893-906.
#   PMID: 16784753   DOI: http://dx.doi.org/10.1016/j.jmb.2006.05.023

# input file:  
#               ../../properties/pdb_prop_geometry.out

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
  
temp = read.table('../../properties/pdb_prop_volume_area.out', header=F, sep="\t", comment.char="", colClasses="character")
pdb_prop_geometry = data.frame( pdb = substr(temp$V6,start=13,stop=18)
                                 , volume = as.numeric(temp$V3)
                                 , area = as.numeric(temp$V5)
                                 , sphericity = 4.8359758620494089221509005399179*(as.numeric(temp$V3)^(2./3.))/as.numeric(temp$V5)
                                 )
pdb_prop_geometry = pdb_prop_geometry[!(pdb_prop_geometry$pdb %in% excluded_pdbs),]
pdb_prop_geometry$pdb  = factor(pdb_prop_geometry$pdb)

temp = read.table('../../properties/res_crd_bf/crd_bf_SC.out', header=T, sep="\t", comment.char="" )
temp         = temp[!(temp$pdb %in% excluded_pdbs),]
temp$pdb     = factor(temp$pdb)

pdb_extent = data.frame()
for ( pdb in levels(pdb_prop_geometry$pdb) )
{
  pdb_crd = as.data.frame( temp[temp$pdb==pdb,] )
  # Now determine the maximum extent of the PDB according to the equation (1) of Lorenz, 1993, " Universality and Cluster Structures in Continuum Models of Percolation "
  max_extent = ( max(pdb_crd$x) - min(pdb_crd$x)
               + max(pdb_crd$y) - min(pdb_crd$y)
               + max(pdb_crd$z) - min(pdb_crd$z)
               ) / 6.
  # Now calculate the Radius of Gyration of the PDB, according to Eqn (2) of Lorenz et al. 1993
  # Or more simply, according to Eqn (45a) Chapter 3.2 of Staufer, Introduction to Percolation Theory
  com_x = mean(pdb_crd$x); com_y = mean(pdb_crd$y); com_z = mean(pdb_crd$z) # center-of-mass coordinates
  dist_from_com = (pdb_crd$x - com_x)^2 + (pdb_crd$y - com_y)^2 + (pdb_crd$z - com_z)^2
  rg_res = sqrt( sum(dist_from_com) / length(dist_from_com) )  # this definition uses the sequence length as the size of the cluster
  rg_vol = sqrt( sum(dist_from_com) / pdb_prop_geometry$volume[pdb_prop_geometry$pdb==pdb] )
  new_row = data.frame( pdb = pdb
                      , nres = length(dist_from_com)
                      , max_extent = max_extent
                      , radius_gyration_nres = rg_res
                      , radius_gyration_vol = rg_vol
                      )
  pdb_extent = rbind ( pdb_extent, new_row )
}


install.packages("MethComp")
library("MethComp")

nres_thresh = 200 # The size at which the protein enters the scaling region according to Lorenz et al. 1993

pdf( "../figures/fractal_dim_max_extent_volume.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.9, 0.9, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot ( log10(pdb_extent$max_extent) , log10(pdb_prop_geometry$volume),  type = 'p', pch=19
     , xlab = expression( "Protein Max. Extent: Log ( R"[m]*" [ Angstroms ] )")
     , ylab = expression( 'Protein Volume: Log ( V [ Angstroms'^3*' ] )' )
     )
reg_gv = Deming( log10(pdb_extent$max_extent), log10(pdb_prop_geometry$volume), boot = TRUE )
abline(reg_gv[1:2], col='red',lwd=2)
graphics.off()

pdf( "../figures/fractal_dim_rad_gyration_nres.pdf", width=4.5, height=4, useDingbats=FALSE )
par( mai=c(0.9, 0.9, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
plot ( log10(pdb_extent$radius_gyration_nres) , log10(pdb_extent$nres),  type = 'p', pch=19
       , xlab = expression( "Radius of Gyration: Log ( R"[g]*" [ Angstroms ] )")
       , ylab = expression( 'Protein Size: Log ( N )' )
     )
reg_gv = Deming( log10(pdb_extent$radius_gyration_nres) , log10(pdb_extent$nres), boot = TRUE )
abline(reg_gv[1:2], col='red',lwd=2)
graphics.off()



plot ( log10(pdb_extent$radius_gyration_vol) , log10(pdb_prop_geometry$volume) )
reg_gv = Deming( log10(pdb_extent$radius_gyration_vol), log10(pdb_prop_geometry$volume), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$max_extent[pdb_extent$nres>nres_thresh]) , log10(pdb_prop_geometry$volume[pdb_extent$nres>nres_thresh]) )
reg_gv = Deming( log10(pdb_extent$max_extent[pdb_extent$nres>nres_thresh]), log10(pdb_prop_geometry$volume[pdb_extent$nres>nres_thresh]), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$max_extent) , log10(pdb_extent$nres) )
reg_gv = Deming( log10(pdb_extent$max_extent), log10(pdb_extent$nres), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$max_extent[pdb_extent$nres>nres_thresh]) , log10(pdb_extent$nres[pdb_extent$nres>nres_thresh]) )
reg_gv = Deming( log10(pdb_extent$max_extent[pdb_extent$nres>nres_thresh]), log10(pdb_extent$nres[pdb_extent$nres>nres_thresh]), boot = TRUE )
abline(reg_gv[1:2])



# RG vs max_extent

plot (log10(pdb_extent$radius_gyration_nres) , log10(pdb_extent$nres) )
reg_gv = Deming( log10(pdb_extent$radius_gyration_nres), log10(pdb_extent$nres), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$radius_gyration_nres[pdb_extent$nres>nres_thresh]) , log10(pdb_extent$nres[pdb_extent$nres>nres_thresh]) )
reg_gv = Deming( log10(pdb_extent$radius_gyration_nres[pdb_extent$nres>nres_thresh]) , log10(pdb_extent$nres[pdb_extent$nres>nres_thresh]), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$radius_gyration_nres) , log10(pdb_extent$max_extent) )
reg_gv = Deming( log10(pdb_extent$radius_gyration_nres) , log10(pdb_extent$max_extent), boot = TRUE )
abline(reg_gv[1:2])

plot (log10(pdb_extent$radius_gyration_vol) , log10(pdb_extent$max_extent) )
reg_gv = Deming( log10(pdb_extent$radius_gyration_vol) , log10(pdb_extent$max_extent), boot = TRUE )
abline(reg_gv[1:2])








all_pdb_prop = read.csv('../tables/all_pdb_prop.csv', header=T)
all_pdb_prop = all_pdb_prop[!(all_pdb_prop$pdb %in% excluded_pdbs),]
all_pdb_prop$variable  = factor(all_pdb_prop$variable)
temp = all_pdb_prop[all_pdb_prop$variable=='nres',]
pdb_prop_geometry = cbind (pdb_prop_geometry, nres = temp$value)

cor(pdb_prop_geometry$sphericity,pdb_prop_geometry$nres,method='sp')
cor(pdb_prop_geometry$volume,pdb_prop_geometry$nres,method='sp')
cor(pdb_prop_geometry$area,pdb_prop_geometry$nres,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC),pdb_prop_geometry$area,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$nres,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_geometry$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_geometry$area,method='sp')
plot(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_geometry$sphericity)
plot(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),log10(pdb_prop_geometry$volume))

cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_geometry$volume,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_geometry$area,method='sp')

library("ppcor")
pcor.test(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_geometry$sphericity,pdb_prop_geometry$volume,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$area,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_geometry$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_geometry$area,method='sp')

cor( cbind( subset(ERscors,select=-c(pdb)) , pdb_prop_geometry[,c("volume","area","sphericity")] ) )
     
cor(abs(ERscors$r.wcnSC.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.rsa.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.rsa.r4sJC),pdb_prop_geometry$volume,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$area,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_geometry$volume,method='sp')

cor(abs(ERscors$r.ddgent.r4sJC),pdb_prop_geometry$sphericity,method='sp')
cor(abs(ERscors$r.ddgent.r4sJC),pdb_prop_geometry$sphericity,method='sp')



# hist(pdb_prop_geometry$sphericity)
# cor(subset(pdb_prop_geometry,select=c(volume,area,sphericity)),method='sp')
# plot(log10(pdb_prop_geometry$volume),log10(pdb_prop_geometry$area))
# plot(log10(pdb_prop_geometry$volume),pdb_prop_geometry$sphericity)
# plot(log10(pdb_prop_geometry$area),pdb_prop_geometry$sphericity)
  

# The following are some random calculations and thoughts. have to figure out what they are later

#  scaling_data = res_prop_voroSC[res_prop_voroSC$VSCvolume_change_diff == 0,]
#  #plot(1.4*(log10(scaling_data$VSCarea)-log10(4.8359758620494089221509005399179)),log10(scaling_data$VSCvolume))
#  plot(log10(scaling_data$VSCarea),log10(scaling_data$VSCvolume/scaling_data$VSCarea^1.25))
#  cor(log10(scaling_data$VSCarea),log10(scaling_data$VSCvolume/scaling_data$VSCarea^1.25),method='sp')
#  #x = 1.4*(log10(scaling_data$VSCarea)-log10(4.8359758620494089221509005399179))
#  x = log10(scaling_data$VSCarea)
#  summary(lm(log10(scaling_data$VSCvolume)~ x))
#  #plot(lm(log10(scaling_data$VSCvolume)~ x))
#  lx = c(1,40)
#  ly = 0.8807092*lx + 0.3023839
#  lines(lx,ly,col='red')

#H-clustering of correlation matrix:

#  all_pdb_prop_cormat = read.csv("../tables/all_pdb_prop_cormat.csv", header=T )
#  row.names(all_pdb_prop_cormat) = all_pdb_prop_cormat$X
#  all_pdb_prop_cormat = subset(all_pdb_prop_cormat, select = -c(X))
#  pdf( "../figures/cormat_hcluster.pdf", width=36, height=16, useDingbats=FALSE )
#  plot( hclust(dist(abs(all_pdb_prop_cormat)))
#      , xlab = 'Hierarchical Clustering of the Spearman Correlation Matrix')
#  dev.off()
#  
#  all_pdb_prop_cormat = read.csv("../tables/all_pdb_prop_cormat.csv", header=T )
#  row.names(all_pdb_prop_cormat) = all_pdb_prop_cormat$X
#  all_pdb_prop_cormat = subset(all_pdb_prop_cormat, select = -c(X))
#  pdf( "../figures/cormat_squared_hcluster.pdf", width=36, height=16, useDingbats=FALSE )
#  plot( hclust(dist(all_pdb_prop_cormat*all_pdb_prop_cormat))
#        , xlab = 'Hierarchical Clustering of the Spearman Correlation Matrix Squared')
#  dev.off()