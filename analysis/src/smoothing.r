# Amir Shahmoradi, Wednesday 3:23 PM, Sep 24 2014, Wilke Lab, ICMB, UT Austin

# install.packages('zoo')
# library('zoo')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

#volume_smooth = filter(res_prop_all_ordered$volume, ma1000)
#plot(res_prop_all_ordered$wcnSC,volume_smooth, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )
#plot(smoothed$wcnSC_mean,smoothed$vol_mean, type='l')

#smoothed = data.frame(vol_mean = rollmean(res_prop_all_ordered$volume, k=2000, no.pad = FALSE),
#                      wcnSC_mean = rollmean(res_prop_all_ordered$wcnSC, k=2000, no.pad = FALSE)
#                      )

# The following is short version of the most useful residue properties that I have found so far.
res_prop_concise = data.frame(seqent        = res_prop_elj$seqent,
                              ddgent        = res_prop_elj$ddgent,
                              rsa           = res_prop_dssp$rsa,
                              hbe           = res_prop_dssp$hbe_mean,
                              hpskd         = res_prop_hps$hpskd,
                              hpsww         = res_prop_hps$hpsww,
                              hpshh         = res_prop_hps$hpshh,
                              wcnSC         = res_prop_wcn_bf$wcnSC,
                              bfSC          = res_prop_wcn_bf$bfSC,
                              resvol        = res_prop_voroSC$resvol,
                              vnvertices    = res_prop_voroSC$vnvertices,
                              vnedges       = res_prop_voroSC$vnedges,
                              vnfaces       = res_prop_voroSC$vnfaces,
                              vedge         = res_prop_voroSC$vedge_length_total,
                              varea         = res_prop_voroSC$varea,
                              vvolume       = res_prop_voroSC$vvolume,
                              veccentricity = res_prop_voroSC$veccentricity
                              )

cormat = cor(res_prop_concise, method='spearman')
write.csv( cormat, "../tables/res_prop_cormat.csv", row.names=T )

# Now since the three Voronoi quantities vnvortices, vnedges and vnfaces happen to be exactly the same as seen in the cormat calculated above, I am going to remove them from the data set, in addition to resvol which is not as informative.
res_prop_concise = subset(res_prop_concise, select = -c(vnvertices,vnedges,resvol))
# The following is an ordered list, in agreement with the column names of the above data frame.
varnames_long = c('Sequence Entropy' , 'ddG Entropy' , 'Relative Solvent Accessibility' , 'Hydrogen Bond Energy' ,
                  'Hydrophobicity Scale (KD)' , 'Hydrophobicity Scale (WW)' , 'Hydrophobicity Scale (HH)' , 'Side-Chain Contact Number' ,
                  'Average Side-Chain B-factor' , 'Voronoi Cell Faces' , 'Voronoi Cell Edge length' , 'Voronoi Cell Surface Area' ,
                  'Voronoi Cell Volume' , 'Voronoi Cell Eccentricity' )

varnames_short = colnames(res_prop_concise)

for (i in 1:length(varnames_short))

  cat(i)
  res_prop_ordered = res_prop_concise[with(res_prop_concise, order(varnames_short[i])),]

test = rollapply(res_prop_ordered, width = 1000, FUN = mean)
test = data.frame(test)

#test_quantile = rollapply(cbind(res_prop_all_ordered$wcnSC,res_prop_all_ordered$volume, res_prop_all_ordered$rsa), width = 1000, FUN = quantile)
#test = data.frame(test_quantile)
#View(test)

plot(test$V,test$V3, type='l')
  
length(res_prop_all_ordered$wcnSC)
length(vol_mean$vol_mean)
