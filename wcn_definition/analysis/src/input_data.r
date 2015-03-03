# This code takes in data generated from my Fortran codes that I wrote in seach of the best performing definitions of WCN for individual structures compared to the original definitions. This data will be further used by other R codes for subsequent analysis.
# Amir Shahmoradi, Friday 10:22 PM, July 18, 2014, Wilke Lab, iCMB, UT Austin

# install.packages('matrixStats')  # needed for functions colMedians
# library(matrixStats)

#setwd('C:/Users/Amir/Documents/GitHub/cordiv/wcn_definition/analysis/src')   # ATTN: Set this to directory where this R code exists


# INPUT DATA FOR POWER-LAW WCN

    sum_wcnp_bfac   = read.table('../../pwrl/bfactor/sum_wcnp_bfac.out',header=T)
    sum_wcnp_seqent = read.table('../../pwrl/seqent/sum_wcnp_seqent.out',header=T)
    exp_wcnp_bfac   = read.table('../../pwrl/bfactor/exp_wcnp_bfac.out',header=F)
    exp_wcnp_seqent = read.table('../../pwrl/seqent/exp_wcnp_seqent.out',header=F)
    
    wcnp_bfac = data.frame()
    for (i in names(exp_wcnp_bfac[,-1])){
       row = data.frame(parameter          = exp_wcnp_bfac[[i]][1],
                        mean_sp       = mean(abs(exp_wcnp_bfac[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcnp_bfac[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcnp_bfac[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcnp_bfac[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcnp_bfac[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcnp_bfac[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcnp_bfac[[i]][-1]))
                        )
       wcnp_bfac = rbind(wcnp_bfac,row)
       }
    
    rownames(wcnp_bfac) = NULL
    
    wcnp_seqent = data.frame()
    for (i in names(exp_wcnp_seqent[,-1])){
       row = data.frame(parameter          = exp_wcnp_seqent[[i]][1],
                        mean_sp       = mean(abs(exp_wcnp_seqent[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcnp_seqent[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcnp_seqent[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcnp_seqent[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcnp_seqent[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcnp_seqent[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcnp_seqent[[i]][-1]))
                        )
       wcnp_seqent = rbind(wcnp_seqent,row)
       }
    
    rownames(wcnp_seqent) = NULL


# INPUT DATA FOR GAUSSIAN WCN

    sum_wcng_bfac   = read.table('../../gaussian/bfactor/sum_wcng_bfac.out',header=T)
    sum_wcng_seqent = read.table('../../gaussian/seqent/sum_wcng_seqent.out',header=T)
    exp_wcng_bfac   = read.table('../../gaussian/bfactor/exp_wcng_bfac.out',header=F)
    exp_wcng_seqent = read.table('../../gaussian/seqent/exp_wcng_seqent.out',header=F)
    
    wcng_bfac = data.frame()
    for (i in names(exp_wcng_bfac[,-1])){
       row = data.frame(parameter     = exp_wcng_bfac[[i]][1],
                        mean_sp       = mean(abs(exp_wcng_bfac[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcng_bfac[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcng_bfac[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcng_bfac[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcng_bfac[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcng_bfac[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcng_bfac[[i]][-1]))
                        )
       wcng_bfac = rbind(wcng_bfac,row)
       }
    
    rownames(wcng_bfac) = NULL
    
    wcng_seqent = data.frame()
    for (i in names(exp_wcng_seqent[,-1])){
       row = data.frame(parameter     = exp_wcng_seqent[[i]][1],
                        mean_sp       = mean(abs(exp_wcng_seqent[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcng_seqent[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcng_seqent[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcng_seqent[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcng_seqent[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcng_seqent[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcng_seqent[[i]][-1]))
                        )
       wcng_seqent = rbind(wcng_seqent,row)
       }
    
    rownames(wcng_seqent) = NULL


# INPUT DATA FOR EXPONENTIAL WCN

    sum_wcne_bfac   = read.table('../../exponential/bfactor/sum_wcne_bfac.out',header=T)
    sum_wcne_seqent = read.table('../../exponential/seqent/sum_wcne_seqent.out',header=T)
    exp_wcne_bfac   = read.table('../../exponential/bfactor/exp_wcne_bfac.out',header=F)
    exp_wcne_seqent = read.table('../../exponential/seqent/exp_wcne_seqent.out',header=F)
    
    wcne_bfac = data.frame()
    for (i in names(exp_wcne_bfac[,-1])){
       row = data.frame(parameter     = exp_wcne_bfac[[i]][1],
                        mean_sp       = mean(abs(exp_wcne_bfac[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcne_bfac[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcne_bfac[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcne_bfac[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcne_bfac[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcne_bfac[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcne_bfac[[i]][-1]))
                        )
       wcne_bfac = rbind(wcne_bfac,row)
       }
    
    rownames(wcne_bfac) = NULL
    
    wcne_seqent = data.frame()
    for (i in names(exp_wcne_seqent[,-1])){
       row = data.frame(parameter     = exp_wcne_seqent[[i]][1],
                        mean_sp       = mean(abs(exp_wcne_seqent[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcne_seqent[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcne_seqent[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcne_seqent[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcne_seqent[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcne_seqent[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcne_seqent[[i]][-1]))
                        )
       wcne_seqent = rbind(wcne_seqent,row)
       }
    
    rownames(wcne_seqent) = NULL


# INPUT DATA FOR DENSITY CN
    
    sum_wcnd_bfac   = read.table('../../density/bfactor/sum_wcnd_bfac.out',header=T)
    sum_wcnd_seqent = read.table('../../density/seqent/sum_wcnd_seqent.out',header=T)
    exp_wcnd_bfac   = read.table('../../density/bfactor/exp_wcnd_bfac.out',header=F)
    exp_wcnd_seqent = read.table('../../density/seqent/exp_wcnd_seqent.out',header=F)
    
    wcnd_bfac = data.frame()
    for (i in names(exp_wcnd_bfac[,-1])){
       row = data.frame(parameter     = exp_wcnd_bfac[[i]][1],
                        mean_sp       = mean(abs(exp_wcnd_bfac[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcnd_bfac[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcnd_bfac[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcnd_bfac[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcnd_bfac[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcnd_bfac[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcnd_bfac[[i]][-1]))
                        )
       wcnd_bfac = rbind(wcnd_bfac,row)
       }
    
    rownames(wcnd_bfac) = NULL
    
    wcnd_seqent = data.frame()
    for (i in names(exp_wcnd_seqent[,-1])){
       row = data.frame(parameter     = exp_wcnd_seqent[[i]][1],
                        mean_sp       = mean(abs(exp_wcnd_seqent[[i]][-1])),
                        median_sp     = quantile(abs(as.vector(exp_wcnd_seqent[[i]][-1])), probs = 0.50),
                        quantile05_sp = quantile(abs(as.vector(exp_wcnd_seqent[[i]][-1])), probs = 0.05),
                        quantile25_sp = quantile(abs(as.vector(exp_wcnd_seqent[[i]][-1])), probs = 0.25),
                        quantile75_sp = quantile(abs(as.vector(exp_wcnd_seqent[[i]][-1])), probs = 0.75),
                        quantile95_sp = quantile(abs(as.vector(exp_wcnd_seqent[[i]][-1])), probs = 0.95),
                        stdev_sp      = sd(as.vector(exp_wcnd_seqent[[i]][-1]))
                        )
       wcnd_seqent = rbind(wcnd_seqent,row)
       }
    
    rownames(wcnd_seqent) = NULL





###  plot(test$parameter,abs(test$mean_sp),type='l', col='black', ylim=c(0.,0.9))
###  lines(test$parameter,abs(test$median_sp),col='red')
###  lines(test$parameter,abs(test$quantile05_sp),col='green')
###  lines(test$parameter,abs(test$quantile25_sp),col='green')
###  lines(test$parameter,abs(test$quantile75_sp),col='green')
###  lines(test$parameter,abs(test$quantile95_sp),col='green')

# test = data.frame(quantile5  = quantile(as.matrix(exp_wcnp_bfac[-1,-1])[,1], probs = 0.05)) #['5%'],
#                   median     = quantile(as.matrix(exp_wcnp_bfac[-1,-1])[,1], probs = c(0.05, 0.5, 0.95))['50%'],
#                   quantile95 = quantile(as.matrix(exp_wcnp_bfac[-1,-1])[,1], probs = c(0.05, 0.5, 0.95))['95%']
#                   )
# 
#wcnp            = data.frame(parameter      = as.vector(t(exp_wcnp_bfac[1,-1])),
#                             mean_sp   = colMeans(exp_wcnp_bfac[-1,-1]),
#                             median_sp = colMedians(as.matrix(exp_wcnp_bfac[-1,-1]))
#                             )