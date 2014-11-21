# This file reads all data for the ASAP project, in order to compare the correlation strengths of seqent-StructuralVariables with the variance of seqent for each pdb.
# Amir Shahmoradi, Thursday 2:08 PM, Oct 1, 2014, ICMB, UT Austin

#INPUT FILES:

setwd('C:/Users/Amir/Documents/GitHub/structural_prediction_of_ER/')

data_1RD8_AB = read.csv( "correlation_analysis/combined_data/data_1RD8_AB.csv" )
data_2FP7_B  = read.csv( "correlation_analysis/combined_data/data_2FP7_B.csv" )
data_2Z83_A  = read.csv( "correlation_analysis/combined_data/data_2Z83_A.csv" )
data_2JLY_A  = read.csv( "correlation_analysis/combined_data/data_2JLY_A.csv" )
data_3GOL_A  = read.csv( "correlation_analysis/combined_data/data_3GOL_A.csv" )
data_3LYF_A  = read.csv( "correlation_analysis/combined_data/data_3LYF_A.csv" ) 
data_4AQF_B  = read.csv( "correlation_analysis/combined_data/data_4AQF_B.csv" )
data_4GHA_A  = read.csv( "correlation_analysis/combined_data/data_4GHA_A.csv" )
data_4IRY_A  = read.csv( "correlation_analysis/combined_data/data_4IRY_A.csv" )
#data_3GSZ_A  = read.csv( "correlation_analysis/combined_data/data_3GSZ_A.csv" )
#data_3I5K_A  = read.csv( "correlation_analysis/combined_data/data_3I5K_A.csv" )

ASAP_res_prop = rbind(data_1RD8_AB, data_2FP7_B, data_2Z83_A, data_2JLY_A, data_3GOL_A, data_3LYF_A, data_4AQF_B, data_4GHA_A, data_4IRY_A)
ASAP_res_prop$protein = factor(ASAP_res_prop$protein)
ASAP_pdb_prop = data.frame()
counter = 0
for (pdb in levels(ASAP_res_prop$protein))
{
  counter = counter + 1
  cat (counter[[1]][1], ' ', pdb, '\n')
  pdb_data = ASAP_res_prop[ASAP_res_prop$protein == pdb,]
  x = cor.test(pdb_data$entropy, pdb_data$rsa_avg_md, method = 'spearman')
  r.seqent.rsa = x$estimate
  x = cor.test(pdb_data$entropy, pdb_data$wcn_avg_md, method = 'spearman')
  r.seqent.wcnCA = x$estimate
  x = cor.test(pdb_data$entropy, pdb_data$bfca, method = 'spearman')
  r.seqent.bfCA = x$estimate
  row = data.frame( pdb  = pdb,
                    nres = length(pdb_data$entropy),
                    r.seqent.rsa   = r.seqent.rsa,
                    r.seqent.wcnCA = r.seqent.wcnCA,
                    r.seqent.bfCA  = r.seqent.bfCA,
                    sd.seqent      = sd(pdb_data$entropy),
                    mean.seqent    = mean(pdb_data$entropy),
                    median.seqent  = median(pdb_data$entropy),
                    sd.wcnCA       = sd(pdb_data$wcn_avg_md),
                    mean.wcnCA     = mean(pdb_data$wcn_avg_md),
                    median.wcnCA   = median(pdb_data$wcn_avg_md),
                    sd.rsa         = sd(pdb_data$rsa_avg_md),
                    mean.rsa       = mean(pdb_data$rsa_avg_md),
                    median.rsa     = median(pdb_data$rsa_avg_md)
                    )
  
  ASAP_pdb_prop = rbind(ASAP_pdb_prop,row)
}

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

write.csv(ASAP_pdb_prop, "../tables/ASAP_pdb_prop.csv", row.names = F)


