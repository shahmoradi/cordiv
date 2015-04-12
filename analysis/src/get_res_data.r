# Last updated by Amir Shahmoradi, Monday 1:51 PM, Feb 2 2015, Wilke Lab, ICMB, UT Austin

# This R script first reads in all residue properties for the entire dataset of 209 proteins from different files into R dataframes.
# input files:  
#               ../../elj_pdb_entropies.in
#               ../../jec_pdb_r4s.csv
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out
#               ../../properties/res_prop_voronoiSC.out
#               ../../properties/res_prop_voronoiAA.out
#               ../../properties/res_prop_voronoiCA.out
#               ../../properties/res_prop_voronoiCB.out
#               ../../properties/res_prop_voronoiN.out
#               ../../properties/res_prop_voronoiC.out
#               ../../properties/res_prop_voronoiO.out

# First read all residue-level data from different tables into R:

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

  excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
  
  # site-specific entropies
  res_prop_elj         = read.table('../../elj_pdb_entropies.in', header=T)
  res_prop_elj         = res_prop_elj[!(res_prop_elj$pdb %in% excluded_pdbs),]
  res_prop_elj$pdb     = factor(res_prop_elj$pdb)
  
  # Rate4Site site-specific evolutionary rates
  res_prop_jec         = read.csv('../../jec_pdb_r4s.csv', header=T)
  res_prop_jec         = res_prop_jec[!(res_prop_jec$pdb %in% excluded_pdbs),]
  res_prop_jec$pdb     = factor(res_prop_jec$pdb)
  
  # Rate4Site site-specific evolutionary rates
  res_prop_jec_ddg     = read.csv('../../jec_ddg_209_monomers.csv', header=T)
  res_prop_jec_ddg     = res_prop_jec_ddg[!(res_prop_jec_ddg$pdb %in% excluded_pdbs),]
  res_prop_jec_ddg$pdb = factor(res_prop_jec_ddg$pdb)
  
  # site-specific Hydrophobicity scales
  res_prop_hps         = read.table('../../properties/res_prop_hps.out', header=T)
  res_prop_hps         = res_prop_hps[!(res_prop_hps$pdb %in% excluded_pdbs),]
  res_prop_hps$pdb     = factor(res_prop_hps$pdb)
  
  # site-specific properties from DSSP software
  res_prop_dssp        = read.table('../../properties/res_prop_dssp.out', header=T)
  res_prop_dssp        = res_prop_dssp[!(res_prop_dssp$pdb %in% excluded_pdbs),]
  res_prop_dssp$pdb    = factor(res_prop_dssp$pdb)
  
  # WCN (power_law kernel) and B factors for different atomic representations of each site
  res_prop_wcn_bf      = read.table('../../properties/res_prop_wcn_bf.out', header=T)
  res_prop_wcn_bf      = res_prop_wcn_bf[!(res_prop_wcn_bf$pdb %in% excluded_pdbs),]
  res_prop_wcn_bf$pdb  = factor(res_prop_wcn_bf$pdb)
  
  # Voronoi cell properties, assuming geometric center of Side Chains as the set of protein corrdinates
  res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
  res_prop_voroSC      = cbind(res_prop_voroSC, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC$VSCvolume^(2./3.))/res_prop_voroSC$VSCarea)
  res_prop_voroSC$VSCmodified_sphericity = res_prop_voroSC$VSCsphericity
  res_prop_voroSC$VSCmodified_sphericity[res_prop_voroSC$VSCvolume_change_diff != 0] = -res_prop_voroSC$VSCsphericity[res_prop_voroSC$VSCvolume_change_diff != 0]
  res_prop_voroSC = res_prop_voroSC[!(res_prop_voroSC$pdb %in% excluded_pdbs),]
  res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
  
  # Voronoi cell properties, assuming geometric center of Amino Acids as the set of protein corrdinates
  res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
  res_prop_voroAA      = cbind(res_prop_voroAA, VAAsphericity = 4.8359758620494089221509005399179*(res_prop_voroAA$VAAvolume^(2./3.))/res_prop_voroAA$VAAarea)
  res_prop_voroAA$VAAmodified_sphericity = res_prop_voroAA$VAAsphericity
  res_prop_voroAA$VAAmodified_sphericity[res_prop_voroAA$VAAvolume_change_diff != 0] = -res_prop_voroAA$VAAsphericity[res_prop_voroAA$VAAvolume_change_diff != 0]
  res_prop_voroAA = res_prop_voroAA[!(res_prop_voroAA$pdb %in% excluded_pdbs),]
  res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)
  
  # Voronoi cell properties, assuming the backbone CA atoms as the set of protein corrdinates
  res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
  res_prop_voroCA      = cbind(res_prop_voroCA, VCAsphericity = 4.8359758620494089221509005399179*(res_prop_voroCA$VCAvolume^(2./3.))/res_prop_voroCA$VCAarea)
  res_prop_voroCA$VCAmodified_sphericity = res_prop_voroCA$VCAsphericity
  res_prop_voroCA$VCAmodified_sphericity[res_prop_voroCA$VCAvolume_change_diff != 0] = -res_prop_voroCA$VCAsphericity[res_prop_voroCA$VCAvolume_change_diff != 0]
  res_prop_voroCA = res_prop_voroCA[!(res_prop_voroCA$pdb %in% excluded_pdbs),]
  res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)
  
  # Voronoi cell properties, assuming the side chain CB atoms as the set of protein corrdinates
  res_prop_voroCB      = read.table('../../properties/res_prop_voronoiCB.out', header=T)
  res_prop_voroCB      = cbind(res_prop_voroCB, VCBsphericity = 4.8359758620494089221509005399179*(res_prop_voroCB$VCBvolume^(2./3.))/res_prop_voroCB$VCBarea)
  res_prop_voroCB$VCBmodified_sphericity = res_prop_voroCB$VCBsphericity
  res_prop_voroCB$VCBmodified_sphericity[res_prop_voroCB$VCBvolume_change_diff != 0] = -res_prop_voroCB$VCBsphericity[res_prop_voroCB$VCBvolume_change_diff != 0]
  res_prop_voroCB = res_prop_voroCB[!(res_prop_voroCB$pdb %in% excluded_pdbs),]
  res_prop_voroCB$pdb  = factor(res_prop_voroCB$pdb)
  
  # Voronoi cell properties, assuming the backbone N atoms as the set of protein corrdinates
  res_prop_voroN      = read.table('../../properties/res_prop_voronoiN.out', header=T)
  res_prop_voroN      = cbind(res_prop_voroN, VNsphericity = 4.8359758620494089221509005399179*(res_prop_voroN$VNvolume^(2./3.))/res_prop_voroN$VNarea)
  res_prop_voroN$VNmodified_sphericity = res_prop_voroN$VNsphericity
  res_prop_voroN$VNmodified_sphericity[res_prop_voroN$VNvolume_change_diff != 0] = -res_prop_voroN$VNsphericity[res_prop_voroN$VNvolume_change_diff != 0]
  res_prop_voroN = res_prop_voroN[!(res_prop_voroN$pdb %in% excluded_pdbs),]
  res_prop_voroN$pdb  = factor(res_prop_voroN$pdb)
  
  # Voronoi cell properties, assuming the backbone C atoms as the set of protein corrdinates
  res_prop_voroC      = read.table('../../properties/res_prop_voronoiC.out', header=T)
  res_prop_voroC      = cbind(res_prop_voroC, VCsphericity = 4.8359758620494089221509005399179*(res_prop_voroC$VCvolume^(2./3.))/res_prop_voroC$VCarea)
  res_prop_voroC$VCmodified_sphericity = res_prop_voroC$VCsphericity
  res_prop_voroC$VCmodified_sphericity[res_prop_voroC$VCvolume_change_diff != 0] = -res_prop_voroC$VCsphericity[res_prop_voroC$VCvolume_change_diff != 0]
  res_prop_voroC = res_prop_voroC[!(res_prop_voroC$pdb %in% excluded_pdbs),]
  res_prop_voroC$pdb  = factor(res_prop_voroC$pdb)
  
  # Voronoi cell properties, assuming the backbone O atoms as the set of protein corrdinates
  res_prop_voroO      = read.table('../../properties/res_prop_voronoiO.out', header=T)
  res_prop_voroO      = cbind(res_prop_voroO, VOsphericity = 4.8359758620494089221509005399179*(res_prop_voroO$VOvolume^(2./3.))/res_prop_voroO$VOarea)
  res_prop_voroO$VOmodified_sphericity = res_prop_voroO$VOsphericity
  res_prop_voroO$VOmodified_sphericity[res_prop_voroO$VOvolume_change_diff != 0] = -res_prop_voroO$VOsphericity[res_prop_voroO$VOvolume_change_diff != 0]
  res_prop_voroO = res_prop_voroO[!(res_prop_voroO$pdb %in% excluded_pdbs),]
  res_prop_voroO$pdb  = factor(res_prop_voroO$pdb)
  

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