# This R code performs PCR analysis on the set of PDB variables in all_pdb_prop_select_wide.csv.
# Last updated by Amir Shahmoradi, Tuesday 4:12 PM, Aug 28 2014, Wilke Lab, ICMB, UT Austin

#install.packages('pls')
library('pls')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

all_pdb_prop_select_wide = read.csv("../tables/all_pdb_prop_select_wide.csv",header=T)
all_pdb_prop_subset = subset(all_pdb_prop_select_wide, select = -c(pdb,sum.nssb,mean.nssb))  # These two columns are all zero for all proteins and do not carry any valuable infromation.

counter = 0
for (column in colnames(all_pdb_prop_subset))
{
  counter = counter + 1
  cat ('processing column # ', counter[[1]], ' : ', column, '\n')
  pcr_out = pcr(all_pdb_prop_subset[[column]] ~ . , data = all_pdb_prop_subset[,-c(counter[[1]])], y = T)
  
  p = pcr_out$projection       # The projection matrix
  
  PEV = c()  # percentage of explained variance
  
  for ( i in 1:pcr_out$ncomp )
  {
    x <- p[,i]
    p[,i] <- x^2/t(x) %*% x
    PEV = cbind(PEV,100*cor(pcr_out$scores[,i],pcr_out$y)^2)
    #print (i)
  }
  
  rownames(PEV) = 'Explained Variance Percentage'
  colnames(PEV) = colnames(p)
  row.names(p) = row.names(pcr_out$projection)
  
  # The first line of the output is the proportion (fraction) of the explained variance of the regressand. The following lines (rows) in the output each represent loading
  PEV_projection = as.data.frame(rbind(PEV,round(pcr_out$projection,6)))
  PEV_projection_ordered = PEV_projection[, c(order(PEV, decreasing = T))]
  
  # Now arrange the variable loadings from large to small for each PCA component, then substitute the loadings with the corresponding variable name.
  # This way each component column will contain the names of the most contributing variables to that component, instead of the projection matrix elements..
  subset_PEV_projection_ordered_values    = PEV_projection_ordered[-1,]     # remove the PEV line from data, then this will contain the prjection matrix elements, ordered for each column (component seperately). 
  subset_PEV_projection_ordered_variables = PEV_projection_ordered[-1,]     # remove the PEV line from data, the n this will contain the names of the variables in the corresponding order of elements in subset_PEV_projection_ordered_values.
  variable_names = row.names(subset_PEV_projection_ordered_values)
  component_names = colnames(subset_PEV_projection_ordered_values)
  
  pmat_ordered = data.frame(row = 1:pcr_out$ncomp)
  pmat_ordered_colnames = c('row')
  PEV_row = data.frame(row = 1)
  for (i in 1:pcr_out$ncomp)
  {
    loadings_order = order(-abs(subset_PEV_projection_ordered_values[[component_names[i]]]))
    subset_PEV_projection_ordered_values[[component_names[i]]] = subset_PEV_projection_ordered_values[[component_names[i]]][loadings_order]  
    subset_PEV_projection_ordered_variables[[component_names[i]]] = variable_names[loadings_order]
    pmat_ordered = cbind(pmat_ordered, subset_PEV_projection_ordered_variables[[component_names[i]]], subset_PEV_projection_ordered_values[[component_names[i]]])
    cvarname = paste0('comp_',strsplit(component_names[i], split = " ")[[1]][2],'_variable')
    cvalname = paste0('comp_',strsplit(component_names[i], split = " ")[[1]][2],'_value')
    pmat_ordered_colnames = rbind(pmat_ordered_colnames, cvarname, cvalname)
    PEV_row = cbind(PEV_row, 'Explained Variance %' , PEV_projection_ordered[[component_names[i]]][1])
  }
  colnames(pmat_ordered) = pmat_ordered_colnames
  
  #PEV_row = as.data.frame(PEV_row)
  colnames(PEV_row) = colnames(pmat_ordered)
  pmat_ordered = rbind(PEV_row,pmat_ordered)
  pmat_ordered = pmat_ordered[,-1]
  
  filename = paste0('../tables/pcr/',column,'_PEV_projection_ordered.csv')
  write.csv( pmat_ordered , filename, row.names=F )

}