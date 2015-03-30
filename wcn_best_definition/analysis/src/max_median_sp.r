# source('input_data.r')

# WCN-Bfactor

    max_median_bfac = data.frame()
    row = data.frame()
    
    row = cbind(model = 'pwrl', wcnp_bfac[which.max(wcnp_bfac$median_sp),])
    max_median_bfac = rbind(max_median_bfac,row)
    
    row = cbind(model = 'gaussian', wcng_bfac[which.max(wcng_bfac$median_sp),])
    max_median_bfac = rbind(max_median_bfac,row)
    
    row = cbind(model = 'exponential', wcne_bfac[which.max(wcne_bfac$median_sp),])
    max_median_bfac = rbind(max_median_bfac,row)
    
    row = cbind(model = 'density', wcnd_bfac[which.max(wcnd_bfac$median_sp),])
    max_median_bfac = rbind(max_median_bfac,row)


# WCN-Seqent

    max_median_seqent = data.frame()
    row = data.frame()
    
    row = cbind(model = 'pwrl', wcnp_seqent[which.max(wcnp_seqent$median_sp),])
    max_median_seqent = rbind(max_median_seqent,row)
    
    row = cbind(model = 'gaussian', wcng_seqent[which.max(wcng_seqent$median_sp),])
    max_median_seqent = rbind(max_median_seqent,row)
    
    row = cbind(model = 'exponential', wcne_seqent[which.max(wcne_seqent$median_sp),])
    max_median_seqent = rbind(max_median_seqent,row)
    
    row = cbind(model = 'density', wcnd_seqent[which.max(wcnd_seqent$median_sp),])
    max_median_seqent = rbind(max_median_seqent,row)