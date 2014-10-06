# Amir Shahmoradi, Wednesday 3:23 PM, Sep 24 2014, Wilke Lab, ICMB, UT Austin

# install.packages('zoo')
library('zoo')

setwd('C:/Users/Amir/Documents/GitHub/cordiv/analysis/src')

# test whether sphericity correlation with wcnSC and other variables could be improved by replacing it with volume_change for open voronoi cells.

modified_sphericity_data = cbind( res_prop_voroSC, modified_sphericity = res_prop_voroSC$VSCsphericity )
modified_sphericity_data$modified_sphericity[modified_sphericity_data$VSCvolume_change_diff != 0] = -modified_sphericity_data$VSCsphericity[modified_sphericity_data$VSCvolume_change_diff != 0]


cor.test(res_prop_voroSC$VSCvolume_change_diff,
         res_prop_voroSC$VSCsphericity,
         method='sp')
cor.test(res_prop_voroSC$VSCvolume_change_diff[res_prop_voroSC$VSCvolume_change_diff != 0],
         res_prop_voroSC$VSCsphericity[res_prop_voroSC$VSCvolume_change_diff != 0],
         method='sp')

