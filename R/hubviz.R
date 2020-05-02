# main function of PALMER


hubviz <- function(dataset, nsample=nrow(dataset), nitem=ncol(dataset), ndim = 2, niter = 30000, nburn = 5000, nthin = 5, nprint = 100,
                    jump_theta = 1.0, jump_w = 0.2, pr_mean_theta = 0.0, pr_sd_theta = 10.0, pr_mean_w = 0.0, 
                    prior_a = 0.001, prior_b = 0.001, w.ci=0.5, option = FALSE, verbose=TRUE, cores = 1){

  message("\n")
  message("---------------------\n")
  message("Running MCMC... \n")
  message("---------------------\n")
  
  if((niter - nburn) %% nthin == 0){

        output = hubviz_cpp(dataset, nsample, nitem, ndim, niter, nburn, nthin, nprint,
                      jump_theta, jump_w, pr_mean_theta, pr_sd_theta, pr_mean_w, 
                      prior_a, prior_b, option, cores)
    
    nmcmc = as.integer((niter - nburn) / nthin)
    max.address = which.max(output$posterior)
    w.star = output$w[max.address,,]
    w.proc = array(0,dim=c(nmcmc,nitem,ndim))
    
    for(iter in 1:nmcmc){
      w.iter = output$w[iter,,]
      if(iter != max.address) w.proc[iter,,] = procrustes(w.iter,w.star)$X.new
      else w.proc[iter,,] = w.iter
    }
    
    theta.est = apply(output$theta,2,mean)
    
    w.est = matrix(NA,nitem,ndim)
    for(i in 1:nitem){
      for(j in 1:ndim){
        w.est[i,j] = mean(w.proc[,i,j])
      }
    }
    sigma.w = mean(output$sigma_w)
    
    w.hpd = array(NA, dim=c(2,nitem,ndim))
    for(i in 1:nitem){
      for(j in 1:ndim){
        w.mcmc = mcmc(w.proc[,i,j])
        w.hpd[,i,j] = HPDinterval(w.mcmc,prob=w.ci)
      }
    }
    
    result = list(theta=output$theta, theta.estimate=theta.est,
                           w=w.proc, w.estimate=w.est, 
                           sigma.w = output$sigma_w, sigma.w.estimate = sigma.w,
                           accept_theta=output$accept_theta,
                           accept_w=output$accept_w, w.hpd=w.hpd)

  }
  else{
    print("Error: The total size of MCMC sample is not integer")
    return(-999)
  }

   if (verbose) {
   result <- result
   } else {
   result <- list( theta.estimate=theta.est, w.estimate=w.est, sigma.w.estimate = sigma.w, w.hpd=w.hpd )
   }

	methods::new("hubviz",
	             data = dataset,
	             init = list(nsample=nsample, nitem=nitem, ndim = ndim, niter = niter, nburn = nburn, nthin=nthin),
	             result = result
	             )
  }
  