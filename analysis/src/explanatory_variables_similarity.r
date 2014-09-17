# This R script attempts at finding the degree of similarity between the correlations of each of the variables with other viarlabes.
# For example, consider r.seqent.rsa, r.seqent.wcnSC. Then this script, first orders chronologically the set of correlating variables with each of these two variables.
# Then compares the Spearman correlations of the two sets with these two vairables with each other.
# Note that only the orser of the sets are important and therefore, only Spearman or Kendall correlations must be used for this purpose.
# A strong correlation would indicate that almost the same sets of varibales in the same order of significance, influence the diversity of the two correlation: r.seqent.rsa, r.seqent.wcnSC
# Amir Shahmoradi, Wednesday 6:36 PM, Sep 3 2014, iCMB, UT Austin

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# source('get_pdb_data.r') ATTN: don't do this! It will take a few hours to get the data summarized at the pdb level. It has been already sourced and the summary exists in the following file.

all_spcor = read.csv( file = '../tables/correlations/all_correlations.csv', header=T )
all_spcor$var1 = factor(all_spcor$var1)

counter1 = 0

for (var11 in levels(all_spcor$var1))
{
  counter2 = 0
  counter1 = counter1 + 1
  factored_spcor1 = all_spcor[all_spcor$var1 == var11,]
  factored_spcor1 = factored_spcor1[with(factored_spcor1, order(var2)),]
  similarity_spcor = data.frame()
  
  if (length(which(!is.na(factored_spcor1$abs.spearman.r))) != 0)
  {
    for (var12 in levels(all_spcor$var1))
    {
      counter2 = counter2 + 1
      cat (counter1,counter2,var11,var12,'\n')
      
      factored_spcor2 = all_spcor[all_spcor$var1 == var12,]
      factored_spcor2 = factored_spcor2[with(factored_spcor2, order(var2)),]
      
      if (var12 != var11 & length(which(!is.na(factored_spcor2$abs.spearman.r))) != 0)
      {
        x = cor.test( factored_spcor1$abs.spearman.r, factored_spcor2$abs.spearman.r, method='spearman', na.action="na.omit" )
        r = x$estimate
        p = x$p.value
        
        row = data.frame( var1 = var11,
                          var2 = var12,
                          abs.spearman.r = abs(r),
                          spearman.r = r,
                          spearman.p = p
                        )
        
        similarity_spcor = rbind(similarity_spcor,row)
      }
    }
    
    similarity_spcor_ordered = similarity_spcor[with(similarity_spcor, order(-abs.spearman.r)),]
    
    filename = paste0( '../tables/correlations/similarity_measures/',var11,'_similarity.csv' )
    write.csv( similarity_spcor_ordered, file = filename, row.names=F )
  }
}
