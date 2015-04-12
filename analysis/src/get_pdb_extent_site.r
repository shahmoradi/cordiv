# Last updated by Amir Shahmoradi, Saturday 10:29 PM, April 11 2015, Wilke Lab, ICMB, UT Austin

# This R script reads in the average side-chain coordinates of all pdbs in the dataset that are found in ../../properties/res_crd_bf/crd_bf_SC.out

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.

temp = read.table('../../properties/res_crd_bf/crd_bf_SC.out', header=T, sep="\t", comment.char="", colClasses="character")
pdb_prop_volume_area = data.frame( pdb = substr(temp$V6,start=13,stop=18)
                                 , volume = as.numeric(temp$V3)
                                 , area = as.numeric(temp$V5)
                                 , sphericity = 4.8359758620494089221509005399179*(as.numeric(temp$V3)^(2./3.))/as.numeric(temp$V5)
                                 )
pdb_prop_volume_area = pdb_prop_volume_area[!(pdb_prop_volume_area$pdb %in% excluded_pdbs),]
pdb_prop_volume_area$pdb  = factor(pdb_prop_volume_area$pdb)

all_pdb_prop = read.csv('../tables/all_pdb_prop.csv', header=T)
all_pdb_prop = all_pdb_prop[!(all_pdb_prop$pdb %in% excluded_pdbs),]
all_pdb_prop$variable  = factor(all_pdb_prop$variable)
temp = all_pdb_prop[all_pdb_prop$variable=='nres',]
pdb_prop_volume_area = cbind (pdb_prop_volume_area, nres = temp$value)

cor(pdb_prop_volume_area$sphericity,pdb_prop_volume_area$nres,method='sp')
cor(pdb_prop_volume_area$volume,pdb_prop_volume_area$nres,method='sp')
cor(pdb_prop_volume_area$area,pdb_prop_volume_area$nres,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC),pdb_prop_volume_area$area,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$nres,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_volume_area$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_volume_area$area,method='sp')
plot(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),pdb_prop_volume_area$sphericity)
plot(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.wcnCA.r4sJC),log10(pdb_prop_volume_area$volume))

cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_volume_area$volume,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_volume_area$area,method='sp')

library("ppcor")
pcor.test(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaCA.r4sJC),pdb_prop_volume_area$sphericity,pdb_prop_volume_area$volume,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$area,method='sp')

cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_volume_area$volume,method='sp')
cor(abs(ERscors$r.wcnSC.r4sJC)-abs(ERscors$r.rsa.r4sJC),pdb_prop_volume_area$area,method='sp')

cor( cbind( subset(ERscors,select=-c(pdb)) , pdb_prop_volume_area[,c("volume","area","sphericity")] ) )
     
cor(abs(ERscors$r.wcnSC.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.rsa.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.rsa.r4sJC),pdb_prop_volume_area$volume,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$area,method='sp')
cor(abs(ERscors$r.vareaSC.r4sJC),pdb_prop_volume_area$volume,method='sp')

cor(abs(ERscors$r.ddgent.r4sJC),pdb_prop_volume_area$sphericity,method='sp')
cor(abs(ERscors$r.ddgent.r4sJC),pdb_prop_volume_area$sphericity,method='sp')

# hist(pdb_prop_volume_area$sphericity)
# cor(subset(pdb_prop_volume_area,select=c(volume,area,sphericity)),method='sp')
# plot(log10(pdb_prop_volume_area$volume),log10(pdb_prop_volume_area$area))
# plot(log10(pdb_prop_volume_area$volume),pdb_prop_volume_area$sphericity)
# plot(log10(pdb_prop_volume_area$area),pdb_prop_volume_area$sphericity)
  

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