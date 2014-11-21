# I am writing this R script in an attempt to find potential structural differences between viral (ASAP) and cellular (Echave) proteins.
# Last updated by Amir Shahmoradi, Wednesday 12:23 PM, Nov 12 2014, Wilke Lab, ICMB, UT Austin

# input files:  
#               ../../elj_pdb_entropies.in
#               ../../properties/res_prop_HPS_asap.out
#               ../../properties/res_prop_dssp_asap.out
#               ../../properties/res_prop_wcn_bf_asap.out


#install.packages("reshape2")
#library("reshape2")
#library('corrplot')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select = read.csv( "../tables/all_pdb_prop_select.csv", header=T )
all_pdb_prop_select$variable = factor(all_pdb_prop_select$variable)
all_pdb_prop_select_asap = read.csv( "../tables/all_pdb_prop_select_asap.csv", header=T )
all_pdb_prop_select_asap$variable = factor(all_pdb_prop_select_asap$variable)

cellular_viral_pdb_differences = data.frame()

for (vvar in levels(all_pdb_prop_select_asap$variable))
{
  if (vvar %in% levels(all_pdb_prop_select$variable))
  {
    # Echave data (213 cellular pdbs)
    cdata = all_pdb_prop_select[all_pdb_prop_select$variable == vvar,]
    x = quantile(cdata$value, probs = c(0,0.25,0.5,0.75, 1.))
    crow = data.frame(c.var      = vvar,
                      c.mean     = mean(cdata$value),
                      c.median   = x[3],
                      c.sd       = sd(cdata$value),
                      c.min      = x[1],
                      c.q1       = x[2],
                      c.q3       = x[4],
                      c.max      = x[5]
                      )
    
    # Viral data (9 pdbs)
    vdata = all_pdb_prop_select_asap[all_pdb_prop_select_asap$variable == vvar,]
    x = quantile(vdata$value, probs = c(0,0.25,0.5,0.75, 1.))
    vrow = data.frame(v.var      = vvar,
                      v.mean     = mean(vdata$value),
                      v.median   = x[3],
                      v.sd       = sd(vdata$value),
                      v.min      = x[1],
                      v.q1       = x[2],
                      v.q3       = x[4],
                      v.max      = x[5]
                      )
                      
    # difference row
    drow = data.frame(c.v.d.mean   = crow$c.mean - vrow$v.mean,
                      c.v.d.median = crow$c.median - vrow$v.median,
                      c.v.d.signif = abs(crow$c.mean - vrow$v.mean)/crow$c.sd
                      )
    row = cbind(crow,vrow,drow)
    cellular_viral_pdb_differences = rbind(cellular_viral_pdb_differences,row)
  }
}

row.names(cellular_viral_pdb_differences) = c()
cellular_viral_pdb_differences = cellular_viral_pdb_differences[with(cellular_viral_pdb_differences, order(-c.v.d.signif)),]
write.csv(cellular_viral_pdb_differences, "../tables/cellular_viral_pdb_differences.csv", row.names=F)
