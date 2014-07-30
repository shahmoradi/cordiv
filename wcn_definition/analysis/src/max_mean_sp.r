# source('input_data.r')

# WCN-Bfactor

    max_mean_bfac = data.frame()
    row = data.frame()
    
    row = cbind(model = 'pwrl', wcnp_bfac[which.max(wcnp_bfac$mean_sp),])
    max_mean_bfac = rbind(max_mean_bfac,row)
    
    row = cbind(model = 'gaussian', wcng_bfac[which.max(wcng_bfac$mean_sp),])
    max_mean_bfac = rbind(max_mean_bfac,row)
    
    row = cbind(model = 'exponential', wcne_bfac[which.max(wcne_bfac$mean_sp),])
    max_mean_bfac = rbind(max_mean_bfac,row)
    
    row = cbind(model = 'density', wcnd_bfac[which.max(wcnd_bfac$mean_sp),])
    max_mean_bfac = rbind(max_mean_bfac,row)


# WCN-Seqent

    max_mean_seqent = data.frame()
    row = data.frame()
    
    row = cbind(model = 'pwrl', wcnp_seqent[which.max(wcnp_seqent$mean_sp),])
    max_mean_seqent = rbind(max_mean_seqent,row)
    
    row = cbind(model = 'gaussian', wcng_seqent[which.max(wcng_seqent$mean_sp),])
    max_mean_seqent = rbind(max_mean_seqent,row)
    
    row = cbind(model = 'exponential', wcne_seqent[which.max(wcne_seqent$mean_sp),])
    max_mean_seqent = rbind(max_mean_seqent,row)
    
    row = cbind(model = 'density', wcnd_seqent[which.max(wcnd_seqent$mean_sp),])
    max_mean_seqent = rbind(max_mean_seqent,row)