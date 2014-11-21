# This code searches for the similarity of the of variables that describe the correlations coefficients.
# This is complicated, ask Amir to explain for you what this code does.
# Amir Shahmoradi, Sunday 1:40 PM, November 2 2014, Wilke Lab, ICMB, UT Austin

#install.packages('corrplot')
library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')



# first Structure_r4sJC correlations:

all_pdb_prop_cormat_SS_scors = read.csv('../tables/cormat_all_pdb_prop_cormat_Sr4s_scors.csv', header=T)

pdf( "../figures/modulators_similarity_Sr4sJC_scors.pdf", width=4.5, height=4, useDingbats=FALSE )
cormat_all_pdb_prop_cormat_SS_scors = cor(subset(all_pdb_prop_cormat_SS_scors, select = -c(variable,vcompactness,vvolume)), method = 'spearman')
corrplot.mixed(cormat_all_pdb_prop_cormat_SS_scors, lower = "ellipse", upper = "number", order = 'alphabet')
dev.off()

# Now Structure_seqent correlations:

all_pdb_prop_cormat_SS_scors = read.csv('../tables/cormat_all_pdb_prop_cormat_Sseqent_scors.csv', header=T)

pdf( "../figures/modulators_similarity_Sseqent_scors.pdf", width=4.5, height=4, useDingbats=FALSE )
cormat_all_pdb_prop_cormat_SS_scors = cor(subset(all_pdb_prop_cormat_SS_scors, select = -c(variable,vcompactness,vvolume)), method = 'spearman')
corrplot.mixed(cormat_all_pdb_prop_cormat_SS_scors, lower = "ellipse", upper = "number", order = 'alphabet')
dev.off()


# Now create the cor matrix of the main modulators of correlations of wcnSC and varea with r4sJC and seqent:

all_pdb_prop_select_wide_rank = read.csv('../tables/all_pdb_prop_select_wide_rank.csv',header=T)
# selection = subset( all_pdb_prop_select_wide_rank
#                   , select = c(r.r4sJC.wcnSC,r.r4sJC.varea,r.seqent.wcnSC,r.seqent.varea
#                               ,r.r4sJC.seqent,sd.seqent,sd.r4sJC
#                               ,sd.hbe,mean.helix,mean.betas))

selection = data.frame( r.r4s.wcn      = all_pdb_prop_select_wide_rank$r.r4sJC.wcnSC
                      , r.r4s.varea    = all_pdb_prop_select_wide_rank$r.r4sJC.varea
                      , r.se.wcn       = all_pdb_prop_select_wide_rank$r.seqent.wcnSC
                      , r.se.varea     = all_pdb_prop_select_wide_rank$r.seqent.varea
                      , r.r4s.se       = all_pdb_prop_select_wide_rank$r.r4sJC.seqent
                      , sd.se          = all_pdb_prop_select_wide_rank$sd.seqent
                      , sd.r4s         = all_pdb_prop_select_wide_rank$sd.r4sJC
                      , sd.hbe         = all_pdb_prop_select_wide_rank$sd.hbe
                      , mn.helix       = all_pdb_prop_select_wide_rank$mean.helix
                      , mn.betas       = all_pdb_prop_select_wide_rank$mean.betas
                      )

pdf( "../figures/main_modulators_cormat.pdf", width=9, height=8, useDingbats=FALSE )
cormat_select = cor(subset(selection), method = 'spearman')
corrplot.mixed(cormat_select, lower = "ellipse", upper = "number"
               #, tl.pos = 'lt'
               , tl.pos = 'd'
               , tl.cex = 0.95
               #, tl.srt = 45
               #, order = 'hclust'
               #, addrect = 3
               )
dev.off()

# Now Generate a more consice version of the above correlation matrix to be included in the main text of the article only for the relation r4sJC-wcnSC

selection = data.frame( r.r4s.wcn      = all_pdb_prop_select_wide_rank$r.r4sJC.wcnSC
                        #, r.r4s.varea    = all_pdb_prop_select_wide_rank$r.r4sJC.varea
                        #, r.se.wcn       = all_pdb_prop_select_wide_rank$r.seqent.wcnSC
                        #, r.se.varea     = all_pdb_prop_select_wide_rank$r.seqent.varea
                        , r.r4s.se       = all_pdb_prop_select_wide_rank$r.r4sJC.seqent
                        , sd.se          = all_pdb_prop_select_wide_rank$sd.seqent
                        , sd.r4s         = all_pdb_prop_select_wide_rank$sd.r4sJC
                        , sd.hbe         = all_pdb_prop_select_wide_rank$sd.hbe
                        , mn.helix       = all_pdb_prop_select_wide_rank$mean.helix
                        , mn.betas       = all_pdb_prop_select_wide_rank$mean.betas
                        )

pdf( "../figures/main_modulators_cormat_r4s_wcn.pdf", width=6, height=5, useDingbats=FALSE )
cormat_select = cor(subset(selection), method = 'spearman')
source ('mycpm.r')
mycpm(cormat_select, lower = "ellipse", upper = "number"
               #, tl.pos = 'lt'
               , tl.pos = 'd'
               , tl.cex = 0.95
               #, tl.srt = 45
               #, order = 'hclust'
               #, addrect = 3
)
dev.off()


# Now Generate the same correlation matrix to be included in the supplementary matterials of the article only for the relation seqent-wcnSC

selection = data.frame( r.r4s.wcn      = all_pdb_prop_select_wide_rank$r.seqent.wcnSC
                        #, r.r4s.varea    = all_pdb_prop_select_wide_rank$r.r4sJC.varea
                        #, r.se.wcn       = all_pdb_prop_select_wide_rank$r.seqent.wcnSC
                        #, r.se.varea     = all_pdb_prop_select_wide_rank$r.seqent.varea
                        , r.r4s.se       = all_pdb_prop_select_wide_rank$r.r4sJC.seqent
                        , sd.se          = all_pdb_prop_select_wide_rank$sd.seqent
                        , sd.r4s         = all_pdb_prop_select_wide_rank$sd.r4sJC
                        , sd.hbe         = all_pdb_prop_select_wide_rank$sd.hbe
                        , mn.helix       = all_pdb_prop_select_wide_rank$mean.helix
                        , mn.betas       = all_pdb_prop_select_wide_rank$mean.betas
                        )

pdf( "../figures/main_modulators_cormat_seqent_wcn.pdf", width=6, height=5, useDingbats=FALSE )
cormat_select = cor(subset(selection), method = 'spearman')
source ('mycpm.r')
mycpm(cormat_select, lower = "ellipse", upper = "number"
      #, tl.pos = 'lt'
      , tl.pos = 'd'
      , tl.cex = 0.95
      #, tl.srt = 45
      #, order = 'hclust'
      #, addrect = 3
)
dev.off()


# Calculate how much of variance is explained by these variables:

#lfit = lm( r.r4sJC.wcnSC ~ sd.seqent + sd.hbe, data = all_pdb_prop_select_wide  )
#lfit = lm( r.seqent.wcnSC ~ sd.seqent + sd.hbe, data = all_pdb_prop_select_wide )
#summary(lfit)

# corrplot(cormat_all_pdb_prop_cormat_SS_scors, order = "hclust")
# cormat_all_pdb_prop_cormat_SS_scors = cor(abs(subset(all_pdb_prop_cormat_SS_scors, select = -c(variable))), method = 'spearman')
# corrplot(cormat_all_pdb_prop_cormat_SS_scors, order = "FPC")
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.r4sJC.rsa)
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.seqent.wcnSC)
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.bfSC.r4sJC)

