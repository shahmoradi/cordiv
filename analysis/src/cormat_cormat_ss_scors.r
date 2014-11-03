# This code searches for the similarity of the of variables that describe the correlations coefficients.
# This is complicated, ask Amir to explain for you what this code does.
# Amir Shahmoradi, Sunday 1:40 PM, November 2 2014, Wilke Lab, ICMB, UT Austin

install.packages('corrplot')
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


# corrplot(cormat_all_pdb_prop_cormat_SS_scors, order = "hclust")
# cormat_all_pdb_prop_cormat_SS_scors = cor(abs(subset(all_pdb_prop_cormat_SS_scors, select = -c(variable))), method = 'spearman')
# corrplot(cormat_all_pdb_prop_cormat_SS_scors, order = "FPC")
# 
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.r4sJC.rsa)
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.seqent.wcnSC)
# plot(all_pdb_prop_cormat_SS_scors$r.r4sJC.wcnSC,all_pdb_prop_cormat_SS_scors$r.bfSC.r4sJC)

