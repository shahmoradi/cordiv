# Last updated by Amir Shahmoradi, Sunday 6:50 PM, April 12 2015, Wilke Lab, ICMB, UT Austin

# setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

  excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.
  
  # read PDB properties
  pdb_prop_general     = read.csv('../../properties/pdb_prop_general.csv', header=T)
  pdb_prop_general     = data.frame( pdb = paste0(pdb_prop_general$Protein,'_',pdb_prop_general$chain), pdb_prop_general )
  pdb_prop_general     = pdb_prop_general[!(pdb_prop_general$pdb %in% excluded_pdbs),]
  pdb_prop_general$pdb = factor(pdb_prop_general$pdb)
  
  # Now read in B factor values for all residues in all proteins
  crd_bf_AA = read.table('../../properties/res_crd_bf/crd_bf_AA.out', header=T)
  crd_bf_AA = crd_bf_AA[!(crd_bf_AA$pdb %in% excluded_pdbs),]
  crd_bf_AA$pdb = factor(crd_bf_AA$pdb)

  crd_bf_CA = read.table('../../properties/res_crd_bf/crd_bf_CA.out', header=T)
  crd_bf_CA = crd_bf_CA[!(crd_bf_CA$pdb %in% excluded_pdbs),]
  crd_bf_CA$pdb = factor(crd_bf_CA$pdb)

  crd_bf_C = read.table('../../properties/res_crd_bf/crd_bf_C.out', header=T)
  crd_bf_C = crd_bf_C[!(crd_bf_C$pdb %in% excluded_pdbs),]
  crd_bf_C$pdb = factor(crd_bf_C$pdb)
  
  pdb_prop_mean_bfAA = data.frame()
  for ( pdb in levels(pdb_prop_general$pdb) )
  {
    pdb_crd_bf_AA = crd_bf_AA[crd_bf_AA$pdb == pdb,]
    pdb_crd_bf_CA = crd_bf_CA[crd_bf_CA$pdb == pdb,]
    pdb_crd_bf_C  = crd_bf_C[crd_bf_C$pdb == pdb,]
    
    mean.bfAA = mean(pdb_crd_bf_AA$bf)
    mean.bfCA = mean(pdb_crd_bf_CA$bf)
    mean.bfC  = mean(pdb_crd_bf_C$bf)
    mean.std_bfAA = mean(pdb_crd_bf_AA$std_bf)
    mean.difmaxmin_bfAA = mean(pdb_crd_bf_AA$bf_max-pdb_crd_bf_AA$bf_min)
    mean.bfCA_bfAA_ratio = mean(pdb_crd_bf_CA$bf / pdb_crd_bf_AA$bf)
    mean.bfC_bfAA_ratio  = mean(pdb_crd_bf_C$bf / pdb_crd_bf_AA$bf)
    row = data.frame( pdb = pdb
                    , mean.bfAA = mean.bfAA
                    , mean.bfCA = mean.bfCA
                    , mean.bfC  = mean.bfC
                    , mean.std_bfAA = mean.std_bfAA
                    , mean.dif_bfAA = mean.dif_bfAA
                    , mean.bfCA_bfAA_ratio = mean.bfCA_bfAA_ratio
                    , mean.bfC_bfAA_ratio  = mean.bfC_bfAA_ratio
                    )
    pdb_prop_mean_bfAA = rbind(pdb_prop_mean_bfAA, row)
  }
  
  cor(pdb_prop_mean_bfAA$mean.bfAA, pdb_prop_mean_bfAA$mean.bfC, method='sp')
  cor(pdb_prop_mean_bfAA$mean.bfAA, pdb_prop_general$Resolution, method='sp')
  cor(pdb_prop_mean_bfAA$mean.bfC, pdb_prop_general$Resolution, method='sp')
  cor(pdb_prop_mean_bfAA$mean.std_bfAA, pdb_prop_general$Resolution, method='sp')
  cor(pdb_prop_mean_bfAA$mean.bfCA_bfAA_ratio, pdb_prop_general$Resolution, method='sp', use = "na.or.complete")
  cor(pdb_prop_mean_bfAA$mean.bfC_bfAA_ratio, pdb_prop_general$Resolution, method='sp', use = "na.or.complete")
  
  
  pdf( "../figures/mean_bfAA_resolution.pdf", width=4.5, height=4, useDingbats=FALSE )
  par( mai=c(0.9, 0.9, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot ( pdb_prop_general$Resolution , log10(pdb_prop_mean_bfAA$mean.bfAA),  type = 'p', pch=19
         , xlab = "X-ray Crystallography Resolution [ Å ] "
         , ylab = expression( 'Mean Atomic B factor: Log ( BF [ Å'^2*' ] )' )
  )
  graphics.off()
  
  pdf( "../figures/mean_bfC_bfAA_ratio_resol.pdf", width=4.5, height=4, useDingbats=FALSE )
  par( mai=c(0.9, 0.9, 0.05, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot ( pdb_prop_general$Resolution , pdb_prop_mean_bfAA$mean.bfC_bfAA_ratio,  type = 'p', pch=19
         , xlab = "X-ray Crystallography Resolution [ Å ] "
         , ylab = expression( 'B factors Ratio: Mean ( BF'[C]*' / BF'[AA]*' )' )
         , ylim = c(0.8,1.1)
  )
  graphics.off()
  
  reg_gv = Deming( pdb_prop_general$Resolution , pdb_prop_mean_bfAA$mean.bfC_bfAA_ratio, boot = TRUE )
  abline(reg_gv[1:2], col='red',lwd=2)
  
  plot(log10(pdb_prop_mean_bfAA$mean.bfAA), pdb_prop_general$Resolution, pch = 16)
  plot(pdb_prop_general$Resolution, pdb_prop_mean_bfAA$mean.bfAA, pch = 16)
  plot(pdb_prop_general$Resolution, pdb_prop_mean_bfAA$mean.bfCA, pch = 16, ylim=c(0,60) )
  mean1 = mean(pdb_prop_mean_bfAA$mean.bfCA[pdb_prop_general$Resolution < 2.], na.rm = TRUE)
  mean2 = mean(pdb_prop_mean_bfAA$mean.bfCA[pdb_prop_general$Resolution > 2.], na.rm = TRUE)
  slope = mean2 - mean1
  reg_lm <- lm(pdb_prop_mean_bfAA$mean.bfAA~pdb_prop_general$Resolution)
  abline(reg_lm,col='red')
  plot(log10(pdb_prop_mean_bfAA$mean.std_bfAA), pdb_prop_general$Resolution, pch = 16)
  plot(pdb_prop_mean_bfAA$mean.bfCA_bfAA_ratio, pdb_prop_general$Resolution, pch = 16)
  plot(pdb_prop_mean_bfAA$mean.bfC_bfAA_ratio, pdb_prop_general$Resolution, pch = 16)
  
  # check if there is any correlation between resolution and protein size
  cor(pdb_prop_general$Sites, pdb_prop_general$Resolution, method='sp')
