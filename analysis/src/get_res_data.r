# Amir Shahmoradi, Monday 7:48 PM, Oct 13 2014, Wilke Lab, ICMB, UT Austin

# reads all residue properties data from different files into R dataframes.

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

excluded_pdbs = c('1BBS_A','1BS0_A','1DIN_A','1HPL_A')   # These are the 4 PDBs that did not have complete r4s evolutionary rates and are omitted from the dataset to avoid NA values.

res_prop_elj         = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj = res_prop_elj[!(res_prop_elj$pdb %in% excluded_pdbs),]
res_prop_elj$pdb     = factor(res_prop_elj$pdb)

res_prop_jec         = read.csv('../../jec_pdb_r4s.csv', header=T)
res_prop_jec = res_prop_jec[!(res_prop_jec$pdb %in% excluded_pdbs),]
res_prop_jec$pdb     = factor(res_prop_jec$pdb)

res_prop_hps         = read.table('../../properties/res_prop_hps.out', header=T)
res_prop_hps = res_prop_hps[!(res_prop_hps$pdb %in% excluded_pdbs),]
res_prop_hps$pdb     = factor(res_prop_hps$pdb)

res_prop_dssp        = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp = res_prop_dssp[!(res_prop_dssp$pdb %in% excluded_pdbs),]
res_prop_dssp$pdb    = factor(res_prop_dssp$pdb)

res_prop_wcn_bf      = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf = res_prop_wcn_bf[!(res_prop_wcn_bf$pdb %in% excluded_pdbs),]
res_prop_wcn_bf$pdb  = factor(res_prop_wcn_bf$pdb)

res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
res_prop_voroAA      = cbind(res_prop_voroAA, VAAsphericity = 4.8359758620494089221509005399179*(res_prop_voroAA$VAAvolume^(2./3.))/res_prop_voroAA$VAAarea)
res_prop_voroAA$VAAmodified_sphericity = res_prop_voroAA$VAAsphericity
res_prop_voroAA$VAAmodified_sphericity[res_prop_voroAA$VAAvolume_change_diff != 0] = -res_prop_voroAA$VAAsphericity[res_prop_voroAA$VAAvolume_change_diff != 0]
res_prop_voroAA = res_prop_voroAA[!(res_prop_voroAA$pdb %in% excluded_pdbs),]
res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
res_prop_voroCA      = cbind(res_prop_voroCA, VCAsphericity = 4.8359758620494089221509005399179*(res_prop_voroCA$VCAvolume^(2./3.))/res_prop_voroCA$VCAarea)
res_prop_voroCA$VCAmodified_sphericity = res_prop_voroCA$VCAsphericity
res_prop_voroCA$VCAmodified_sphericity[res_prop_voroCA$VCAvolume_change_diff != 0] = -res_prop_voroCA$VCAsphericity[res_prop_voroCA$VCAvolume_change_diff != 0]
res_prop_voroCA = res_prop_voroCA[!(res_prop_voroCA$pdb %in% excluded_pdbs),]
res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC$VSCvolume^(2./3.))/res_prop_voroSC$VSCarea)
res_prop_voroSC$VSCmodified_sphericity = res_prop_voroSC$VSCsphericity
res_prop_voroSC$VSCmodified_sphericity[res_prop_voroSC$VSCvolume_change_diff != 0] = -res_prop_voroSC$VSCsphericity[res_prop_voroSC$VSCvolume_change_diff != 0]
res_prop_voroSC = res_prop_voroSC[!(res_prop_voroSC$pdb %in% excluded_pdbs),]
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)