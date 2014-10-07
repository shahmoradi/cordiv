# This R function takes in all the the residue properties of all proteins from all sorts of analysis, and then selects the  most relevant quantites that will be later used in the analyses. Unlike its ancestor code 'get_res_data.r', this code only selects variables that performed the best among other similar variables definitions. For example, among all Bfactor definitions, only the average over SideChain atoms (bfSC) will be taken and used. This is because apparently, as far as I have searched, residue properties that are based on side chain atoms, correlate best with other residue properties. The only exception to this is the H-bond energy (hbe_mean) of residues which correlates best with CA atom properties (such as bfCA and wcnCA).
# Also, among the three hydrophobicity scales, I am only picking hpshh, which seems to correlate best with other residue properties.

# input files:  
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_hps.out
#               ../../properties/res_prop_dssp.out
#               ../../properties/res_prop_wcn_bf.out
#               ../../properties/res_prop_voroSC.out

# Last updated by Amir Shahmoradi, Thursday 3:40 PM, Aug 28 2014, Wilke Lab, ICMB, UT Austin

#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

res_prop_elj         = read.table('../../elj_pdb_entropies.in', header=T)
res_prop_elj$pdb     = factor(res_prop_elj$pdb)

res_prop_hps         = read.table('../../properties/res_prop_hps.out', header=T)
res_prop_hps$pdb     = factor(res_prop_hps$pdb)

res_prop_dssp        = read.table('../../properties/res_prop_dssp.out', header=T)
res_prop_dssp$pdb    = factor(res_prop_dssp$pdb)

res_prop_wcn_bf      = read.table('../../properties/res_prop_wcn_bf.out', header=T)
res_prop_wcn_bf$pdb  = factor(res_prop_wcn_bf$pdb)

#res_prop_voroAA      = read.table('../../properties/res_prop_voronoiAA.out', header=T)
#res_prop_voroAA$pdb  = factor(res_prop_voroAA$pdb)

#res_prop_voroCA      = read.table('../../properties/res_prop_voronoiCA.out', header=T)
#res_prop_voroCA$pdb  = factor(res_prop_voroCA$pdb)

res_prop_voroSC      = read.table('../../properties/res_prop_voronoiSC.out', header=T)
res_prop_voroSC      = cbind(res_prop_voroSC, VSCsphericity = 4.8359758620494089221509005399179*(res_prop_voroSC$VSCvolume^(2./3.))/res_prop_voroSC$VSCarea)
res_prop_voroSC$VSCmodified_sphericity = res_prop_voroSC$VSCsphericity
res_prop_voroSC$VSCmodified_sphericity[res_prop_voroSC$VSCvolume_change_diff != 0] = -res_prop_voroSC$VSCsphericity[res_prop_voroSC$VSCvolume_change_diff != 0]
res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
# res_prop_voroSC$pdb  = factor(res_prop_voroSC$pdb)
# res_prop_voroSC      = cbind(res_prop_voroSC, VSCmodified_volume = res_prop_voroSC$VSCvolume)
# maxval = max(res_prop_voroSC$VSCvolume)
# res_prop_voroSC$VSCmodified_volume[res_prop_voroSC$VSCvolume_change != 0] = maxval
# res_prop_voroSC$VSCmodified_volume = res_prop_voroSC$VSCmodified_volume + res_prop_voroSC$VSCvolume_change

pdb_prop_from_residue_prop_select = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_elj$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_elj    = res_prop_elj[res_prop_elj$pdb==pdb,] # c('seqent','ddgent')]
  pdb_hps    = res_prop_hps[res_prop_hps$pdb==pdb,] # c('hpskd','hpsww','hpshh')] )
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )
  pdb_wcn_bf = res_prop_wcn_bf[res_prop_wcn_bf$pdb==pdb, ]
  pdb_voroSC = res_prop_voroSC[res_prop_voroSC$pdb==pdb, ]
  
  pdb_temp = data.frame( seqent   = pdb_elj$seqent,
                         ddgent   = pdb_elj$ddgent,
                         hpshh    = pdb_hps$hpshh,
                         rsa      = pdb_dssp$rsa,
                         hbe      = pdb_dssp$hbe_mean,
                         bfSC     = pdb_wcn_bf$bfSC,
                         wcnSC    = pdb_wcn_bf$wcnSC,
                         vnfaces       = pdb_voroSC$VSCnfaces,
                         vedge         = log10(pdb_voroSC$VSCedge_length_total),
                         varea         = log10(pdb_voroSC$VSCarea),
                         vvolume       = log10(pdb_voroSC$VSCvolume),
                         veccentricity = pdb_voroSC$VSCeccentricity,
                         vsphericity   = pdb_voroSC$VSCsphericity,
                         mvsphericity  = pdb_voroSC$VSCmodified_sphericity
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
      row = data.frame(pdb, variable = paste0('sum.',variable1), value = sum(var1$value))       ; pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,row)
      row = data.frame(pdb, variable = paste0('mean.',variable1), value = mean(var1$value))     ; pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,row)
      row = data.frame(pdb, variable = paste0('median.',variable1), value = median(var1$value)) ; pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,row)
      row = data.frame(pdb, variable = paste0('sd.',variable1), value = sd(var1$value))         ; pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,row)
      
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
          pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,row)
        }
      }
  }

}


# Now get Secondary Structure residue data:
pdb_prop_ss = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_dssp$pdb))
{
  counter = counter + 1
  cat( paste(str(counter),pdb) )
  
  pdb_dssp   = res_prop_dssp[res_prop_dssp$pdb==pdb,] # c('asa','rsa','hbe_mean','rss')] )

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

pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,pdb_prop_ss)

# One last step: Since voronoi both sphericity (vsphericity) and modified voronoi sphericity (mvsphericity) are biased quantities by definition, their mean, median, and standard deviations (sd) may not be so meaningful.
# Therefore, in the following lines I am also going to calculate the moments of vsphericity for only those voronoi cells that are closed (vccsphericity), so that the amount bias due to edge effects could be minimized.
pdb_prop_closed_cells = data.frame()    # This dataframe will contain the mean median and variance of sequqence entropy and ddG entropy for each pdb file.
counter = 0

for(pdb in levels(res_prop_voroSC$pdb))
{
  counter = counter + 1
  cat( 'vccsphericity: ', paste(str(counter),pdb), '\n' )
  # Now select only those residues that belong to protein pdb and have closed voronoi cells.
  pdb_voroSC = res_prop_voroSC[(res_prop_voroSC$pdb==pdb & res_prop_voroSC$VSCvolume_change_diff==0),]
  # Calculate potentially important statistical moments of the VSCsphericity:
  row = data.frame(pdb, variable = 'sum.vccsphericity', value = sum(pdb_voroSC$VSCsphericity))       ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'mean.vccsphericity', value = mean(pdb_voroSC$VSCsphericity))     ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'median.vccsphericity', value = median(pdb_voroSC$VSCsphericity)) ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'sd.vccsphericity', value = sd(pdb_voroSC$VSCsphericity))         ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  # Calculate potentially important statistical moments of the VSCnfaces:
  row = data.frame(pdb, variable = 'sum.vccnfaces', value = sum(pdb_voroSC$VSCnfaces))       ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'mean.vccnfaces', value = mean(pdb_voroSC$VSCnfaces))     ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'median.vccnfaces', value = median(pdb_voroSC$VSCnfaces)) ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
  row = data.frame(pdb, variable = 'sd.vccnfaces', value = sd(pdb_voroSC$VSCnfaces))         ; pdb_prop_closed_cells = rbind(pdb_prop_closed_cells,row)
}
pdb_prop_from_residue_prop_select = rbind(pdb_prop_from_residue_prop_select,pdb_prop_closed_cells)


write.csv( pdb_prop_from_residue_prop_select, "../tables/pdb_prop_from_residue_prop_select.csv", row.names=F )

