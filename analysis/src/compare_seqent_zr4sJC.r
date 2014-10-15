# This R code compares the correlation of the two measures of sequence variability (seqent and zr4sJC) with different structural properties.
# Amir Shahmoradi, Tuesday 3:06 PM, Oct 14 2014, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# source('get_res_data.r')

# The following properties will be considered for comparison:
# ddGent , rsa , wcnSC , bfSC , VSCnfaces , VSCedge , VSCarea , VSCvolume , VSCeccentricity , mVSCsphericity , hbe , hpshh

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv", header = T)
seqent_cors = c('r.ddgent.seqent','r.rsa.seqent','r.seqent.wcnSC','r.seqent.vedge','r.seqent.varea','r.seqent.vvolume','r.bfSC.seqent','r.seqent.veccentricity','r.seqent.vsphericity','r.seqent.vnfaces','r.hbe.seqent','r.hpshh.seqent')
r4sJC_cors = c('r.ddgent.r4sJC','r.r4sJC.rsa','r.r4sJC.wcnSC','r.r4sJC.vedge','r.r4sJC.varea','r.r4sJC.vvolume','r.bfSC.r4sJC','r.r4sJC.veccentricity','r.r4sJC.vsphericitym','r.r4sJC.vnfaces','r.hbe.r4sJC','r.hpshh.r4sJC')

x = -1:2.
pdf( '../figures/seqent_vs_zr4sJC.pdf', width=13.5, height=16, useDingbats=FALSE )  
split.screen(c(4,3))

scor_stat = data.frame()

for (i in 1:length(seqent_cors))
{
  cat(i)
  screen(i)
  par( mai=c(0.65, 0.65, 0.1, 0.05), mgp=c(2, 0.5, 0), tck=-0.03 )
  plot(abs(all_pdb_prop_select_wide[[seqent_cors[i]]]),
       abs(all_pdb_prop_select_wide[[r4sJC_cors[i]]]),
       pch = 16,
       xlim = c(0,0.9),
       ylim = c(0,0.9),
       xlab = paste0('ABS ( ',seqent_cors[i],' )'),
       ylab = paste0('ABS ( ',r4sJC_cors[i],' )'),
       cex.lab = 1.4,
       cex.axis = 1.4
       )
  lines(x,x,col='red')
  #lines(x,x+0.1,col='red')
  
  #row = data.frame(variable = )
}
dev.off()

