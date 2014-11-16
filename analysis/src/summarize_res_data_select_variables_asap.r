# This R function takes in the residue properties of all ASAP proteins from DSSP, and PDB files and and combines them into a single dataframe for further analysis at pdb level.

# Last updated by Amir Shahmoradi, Wednesday 3:02 PM, Nov 12 2014, Wilke Lab, ICMB, UT Austin

# input files:  
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps_asap.out
#               ../../properties/res_prop_dssp_asap.out
#               ../../properties/res_prop_wcn_bf_asap.out
#               ../../properties/res_prop_voroSC_asap.out


#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# First get ASAP entropies:

seqent_1RD8_X = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/1RD8_X_HP_entropies.csv')
seqent_1RD8_X = data.frame(pdb = '1RD8_X',seqent_1RD8_X)
seqent_2FP7_B = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/2FP7_B_WNPB_entropies.csv')
seqent_2FP7_B = data.frame(pdb = '2FP7_B',seqent_2FP7_B)
seqent_2JLY_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/2JLY_A_DPH_entropies.csv')
seqent_2JLY_A = data.frame(pdb = '2JLY_A',seqent_2JLY_A)
seqent_2Z83_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/2Z83_A_JEHN_entropies.csv')
seqent_2Z83_A = data.frame(pdb = '2Z83_A',seqent_2Z83_A)
seqent_3GOL_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/3GOL_A_HCP_entropies.csv')
seqent_3GOL_A = data.frame(pdb = '3GOL_A',seqent_3GOL_A)
seqent_3GSZ_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/3GSZ_A_HCP_entropies.csv')
seqent_3GSZ_A = data.frame(pdb = '3GSZ_A',seqent_3GSZ_A)
seqent_3I5K_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/3I5K_A_HCP_entropies.csv')
seqent_3I5K_A = data.frame(pdb = '3I5K_A',seqent_3I5K_A)
seqent_3LYF_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/3LYF_A_RVFVNP_entropies.csv')
seqent_3LYF_A = data.frame(pdb = '3LYF_A',seqent_3LYF_A)
seqent_4AQF_B = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/4AQF_B_CCHFN_entropies.csv')
seqent_4AQF_B = data.frame(pdb = '4AQF_B',seqent_4AQF_B)
seqent_4GHA_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/4GHA_A_MRNABD_entropies.csv')
seqent_4GHA_A = data.frame(pdb = '4GHA_A',seqent_4GHA_A)
seqent_4IRY_A = read.csv('../../rate_calculations/entropy/prelim_entropies/entropies/4IRY_A_INP_entropies.csv')
seqent_4IRY_A = data.frame(pdb = '4IRY_A',seqent_4IRY_A)

res_prop_seqent = rbind( seqent_1RD8_X
                       , seqent_2FP7_B
                       , seqent_2JLY_A
                       , seqent_2Z83_A
                       , seqent_2Z83_A
                       , seqent_3GOL_A
                       , seqent_3GSZ_A
                       , seqent_3I5K_A
                       , seqent_3LYF_A
                       , seqent_4AQF_B
                       , seqent_4GHA_A
                       , seqent_4IRY_A
                       )
res_prop_seqent$pdb = data.frame(res_prop_seqent$pdb)

res_prop_hps_asap        = read.table('../../properties/res_prop_hps_asap.out', header=T)
res_prop_hps_asap$pdb    = factor(res_prop_hps_asap$pdb)
                         
res_prop_dssp_asap       = read.table('../../properties/res_prop_dssp_asap.out', header=T)
res_prop_dssp_asap$pdb   = factor(res_prop_dssp_asap$pdb)

res_prop_wcn_bf_asap     = read.table('../../properties/res_prop_wcn_bf_asap.out', header=T)
res_prop_wcn_bf_asap$pdb = factor(res_prop_wcn_bf_asap$pdb)

#res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA_asap.out', header=T)
#res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

#res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA_asap.out', header=T)
#res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

res_prop_voroSC_asap      = read.table('../../properties/res_prop_voronoiSC_asap.out', header=T)
res_prop_voroSC_asap      = cbind(res_prop_voroSC_asap, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC_asap$VSCvolume^(2./3.))/res_prop_voroSC_asap$VSCarea)
res_prop_voroSC_asap$VSCmodified_sphericity = res_prop_voroSC_asap$VSCsphericity
res_prop_voroSC_asap$VSCmodified_sphericity[res_prop_voroSC_asap$VSCvolume_change_diff != 0] = -res_prop_voroSC_asap$VSCsphericity[res_prop_voroSC_asap$VSCvolume_change_diff != 0]
res_prop_voroSC_asap$pdb  = factor(res_prop_voroSC_asap$pdb)
# res_prop_voroSC_asap$pdb  = factor(res_prop_voroSC_asap$pdb)
# res_prop_voroSC_asap      = cbind(res_prop_voroSC_asap, VSCmodified_volume = res_prop_voroSC_asap$VSCvolume)
# maxval = max(res_prop_voroSC_asap$VSCvolume)
# res_prop_voroSC_asap$VSCmodified_volume[res_prop_voroSC_asap$VSCvolume_change != 0] = maxval
# res_prop_voroSC_asap$VSCmodified_volume = res_prop_voroSC_asap$VSCmodified_volume + res_prop_voroSC_asap$VSCvolume_change

nonviral_pdbs = c('1AJ8_A','1AOR_A','1CTS_A','1MP9_A','3GSZ_A','3I5K_A')   # These are thermophilic proteins, in addition to the 2 viral PDBs that are from the same family as 3GOL_A and therefore redundant.
pdb_prop_from_residue_prop_select_asap = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_dssp_asap$pdb))
{
  if (!(pdb %in% nonviral_pdbs))
  {
    counter = counter + 1
    cat( paste(pdb, str(counter[[1]][1]),'\n') )
    
    pdb_seqent  = res_prop_seqent[res_prop_seqent$pdb==pdb,] # c('seqent','ddgent')]
    #pdb_jec    = res_prop_jec[res_prop_jec$pdb==pdb,] # c('zr4s_JC')]
    pdb_hps    = res_prop_hps_asap[res_prop_hps_asap$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
    pdb_dssp   = res_prop_dssp_asap[res_prop_dssp_asap$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
    pdb_wcn_bf = res_prop_wcn_bf_asap[res_prop_wcn_bf_asap$pdb==pdb, ]
    pdb_voroSC = res_prop_voroSC_asap[res_prop_voroSC_asap$pdb==pdb, ]
    
    pdb_temp = data.frame( seqent   = pdb_seqent$entropy
                         #, ddgent   = pdb_elj$ddgent
                         #, r4sJC    = pdb_jec$r4s_JC
                         #, r4sJCz   = pdb_jec$zr4s_JC
                         , hpshh    = pdb_hps$hpshh
                         , rsa      = pdb_dssp$rsa
                         , hbe      = pdb_dssp$hbe
                         , bfSC     = pdb_wcn_bf$bfSC
                         , wcnSC    = pdb_wcn_bf$wcnSC
                         , vnfaces       = pdb_voroSC$VSCnfaces
                         , vedge         = log10(pdb_voroSC$VSCedge_length_total)
                         , varea         = log10(pdb_voroSC$VSCarea)
                         , vvolume       = log10(pdb_voroSC$VSCvolume)
                         , veccentricity = pdb_voroSC$VSCeccentricity
                         , vsphericity   = pdb_voroSC$VSCsphericity
                         , vsphericitym  = pdb_voroSC$VSCmodified_sphericity
                         )
    
    pdb_long = reshape(pdb_temp, ids = rownames(pdb_temp), varying = colnames(pdb_temp), v.names = 'value', timevar = 'variable', times = colnames(pdb_temp), direction = 'long')
    pdb_long$variable = factor(pdb_long$variable)
    
    counter1 = 0
    
    for (variable1 in levels(pdb_long$variable))
    {
      counter1 = counter1 + 1
      #cat (variable1, '\n')
      var1 = pdb_long[pdb_long$variable == variable1,]
      
      # Calculate potentially important statistical moments of the factored variable:
      if(variable1 != 'r4sJCz')     # There is no information in the moments of the standardized values r4sJCz.
      {
        row = data.frame(pdb, variable = paste0('sum.',variable1), value = sum(var1$value))       ; pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,row)
        row = data.frame(pdb, variable = paste0('mean.',variable1), value = mean(var1$value))     ; pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,row)
        row = data.frame(pdb, variable = paste0('median.',variable1), value = median(var1$value)) ; pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,row)
        row = data.frame(pdb, variable = paste0('sd.',variable1), value = sd(var1$value))         ; pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,row)
      }
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
          pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,row)
        }
      }
    }
  }
}


# Now get Secondary Structure residue data:
pdb_prop_ss = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_dssp_asap$pdb))
{
  if (!(pdb %in% nonviral_pdbs))
  {
    counter = counter + 1
    cat( str(counter[1]), pdb, '\n' )
    
    pdb_dssp   = res_prop_dssp_asap[res_prop_dssp_asap$pdb==pdb,] # c('asa','rsa','hbe','rss')] )
    
    sum.GSS    = length(which(pdb_dssp$rss == 'G'))
    sum.HSS    = length(which(pdb_dssp$rss == 'H'))
    sum.ISS    = length(which(pdb_dssp$rss == 'I'))
    sum.TSS    = length(which(pdb_dssp$rss == 'T'))
    sum.ESS    = length(which(pdb_dssp$rss == 'E'))
    sum.BSS    = length(which(pdb_dssp$rss == 'B'))
    sum.SSS    = length(which(pdb_dssp$rss == 'S'))
    sum.LSS    = length(which(pdb_dssp$rss == 'L'))
               + length(which(pdb_dssp$rss == 'C'))
               + length(which(pdb_dssp$rss == '_'))
    sum.helix  = sum.GSS + sum.HSS + sum.ISS
    sum.betas  = sum.ESS + sum.BSS
    sum.hbdif  = sum.helix + sum.betas
    
    pdb.nres   = length(pdb_dssp$pdb)
    
    mean.GSS   = sum.GSS/pdb.nres
    mean.HSS   = sum.HSS/pdb.nres
    mean.ISS   = sum.ISS/pdb.nres
    mean.TSS   = sum.TSS/pdb.nres
    mean.ESS   = sum.ESS/pdb.nres
    mean.BSS   = sum.BSS/pdb.nres
    mean.SSS   = sum.SSS/pdb.nres
    mean.LSS   = sum.LSS/pdb.nres
    mean.helix = sum.helix/pdb.nres
    mean.betas = sum.betas/pdb.nres
    mean.hbdif = sum.hbdif/pdb.nres
    
    row = data.frame(pdb, variable = 'sum.GSS'   , value = sum.GSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.HSS'   , value = sum.HSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.ISS'   , value = sum.ISS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.TSS'   , value = sum.TSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.ESS'   , value = sum.ESS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.BSS'   , value = sum.BSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.SSS'   , value = sum.SSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.LSS'   , value = sum.LSS   ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.helix' , value = sum.helix ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.betas' , value = sum.betas ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'sum.hbdif' , value = sum.hbdif ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.GSS'  , value = mean.GSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.HSS'  , value = mean.HSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.ISS'  , value = mean.ISS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.TSS'  , value = mean.TSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.ESS'  , value = mean.ESS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.BSS'  , value = mean.BSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.SSS'  , value = mean.SSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.LSS'  , value = mean.LSS  ) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.helix', value = mean.helix) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.betas', value = mean.betas) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
    row = data.frame(pdb, variable = 'mean.hbdif', value = mean.hbdif) ; pdb_prop_ss = rbind(pdb_prop_ss,row)
  }
}

pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,pdb_prop_ss)

# One last step: Since voronoi both sphericity (vsphericity) and modified voronoi sphericity (mvsphericity) are biased quantities by definition, their mean, median, and standard deviations (sd) may not be so meaningful.
# Therefore, in the following lines I am also going to calculate the moments of vsphericity for only those voronoi cells that are closed (vccsphericity), so that the amount bias due to edge effects could be minimized.
pdb_prop_closed_cells = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_voroSC_asap$pdb))
{
  if (!(pdb %in% nonviral_pdbs))
  {
    counter = counter + 1
    cat( 'vccsphericity: ', paste(str(counter),pdb), '\n' )
    # Now select only those residues that belong to protein pdb and have closed voronoi cells.
    pdb_voroSC = res_prop_voroSC_asap[(res_prop_voroSC_asap$pdb==pdb & res_prop_voroSC_asap$VSCvolume_change_diff==0),]
    # Calculate potentially important statistical moments of the VSCsphericity:
    row = data.frame(pdb, variable = 'sum.vsphericitycc', value = sum(pdb_voroSC$VSCsphericity))       ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'mean.vsphericitycc', value = mean(pdb_voroSC$VSCsphericity))     ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'median.vsphericitycc', value = median(pdb_voroSC$VSCsphericity)) ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'sd.vsphericitycc', value = sd(pdb_voroSC$VSCsphericity))         ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    # Calculate potentially important statistical moments of the VSCnfaces:
    row = data.frame(pdb, variable = 'sum.vnfacescc', value = sum(pdb_voroSC$VSCnfaces))       ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'mean.vnfacescc', value = mean(pdb_voroSC$VSCnfaces))     ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'median.vnfacescc', value = median(pdb_voroSC$VSCnfaces)) ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
    row = data.frame(pdb, variable = 'sd.vnfacescc', value = sd(pdb_voroSC$VSCnfaces))         ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  }
}
pdb_prop_from_residue_prop_select_asap = rbind(pdb_prop_from_residue_prop_select_asap,pdb_prop_closed_cells)


write.csv( pdb_prop_from_residue_prop_select_asap, "../tables/pdb_prop_from_residue_prop_select_asap.csv", row.names=F )

