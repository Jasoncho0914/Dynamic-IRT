build_Q = function(obs_sigma_t2, evol_sigma_t2, D = 1){
  
  if(!(D == 1 || D == 2)) stop('build_Q requires D = 1 or D = 2')
  
  T = length(evol_sigma_t2)
  
  # For reference: first and second order difference matrices (not needed below)
  #H1 = bandSparse(T, k = c(0,-1), diag = list(rep(1, T), rep(-1, T)), symmetric = FALSE)
  #H2 = bandSparse(T, k = c(0,-1, -2), diag = list(rep(1, T), c(0, rep(-2, T-1)), rep(1, T)), symmetric = FALSE)
  
  # Quadratic term: can construct directly for D = 1 or D = 2 using [diag(1/obs_sigma_t2, T) + (t(HD)%*%diag(1/evol_sigma_t2, T))%*%HD]
  if(D == 1){
    # D = 1 case:
    Q = Matrix::bandSparse(T, k = c(0,1),
                           diagonals= list(1/obs_sigma_t2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0),
                                           -1/evol_sigma_t2[-1]),
                           symmetric = TRUE)
  } else {
    # D = 2 case:
    Q = Matrix::bandSparse(T, k = c(0,1,2),
                           diagonals= list(1/obs_sigma_t2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0),
                                           c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0))),
                                           1/evol_sigma_t2[-(1:2)]),
                           symmetric = TRUE)
  }
  Q
}
initEvol0 = function(mu0, commonSD = TRUE){
  
  p = length(mu0)
  
  # Common or distinct:
  if(commonSD) {
    sigma_w0 = rep(mean(abs(mu0)), p)
  } else  sigma_w0 = abs(mu0)
  
  # Initialize at 1 for simplicity:
  px_sigma_w0 = rep(1, p)
  
  sigma_00 = px_sigma_00 = 1
  
  list(sigma_w0 = sigma_w0, px_sigma_w0 = px_sigma_w0, sigma_00 = sigma_00, px_sigma_00 = px_sigma_00)
}
sampleBTF = function(y, obs_sigma_t2, evol_sigma_t2, D = 1, loc_obs = NULL, chol0 = NULL){
  
  # Some quick checks:
  if((D < 0) || (D != round(D)))  stop('D must be a positive integer')
  
  if(any(is.na(y))) stop('y cannot contain NAs')
  
  T = length(y)
  Td = T-D
  # Linear term:
  linht = y/obs_sigma_t2
  
  # Quadratic terms and solutions are computed differently, depending on D:
  
  if(D == 0){
    # Special case: no differencing
    
    # Posterior SDs and posterior means:
    postSD = 1/sqrt(1/obs_sigma_t2 + 1/evol_sigma_t2)
    postMean = (linht)*postSD^2
    
    # Sample the states:
    mu = rnorm(n = T, mean = postMean, sd = postSD)
    
  } else {
    
    # New sampler, based on spam package:
    QHt_Matrix = build_Q(obs_sigma_t2 = obs_sigma_t2, evol_sigma_t2 = evol_sigma_t2, D = D)
    
    if(!is.null(chol0)){
      
      # Sample the states:
      mu = matrix(rmvnorm.canonical(n = 1,
                                    b = linht,
                                    Q = as.spam.dgCMatrix(as(QHt_Matrix, "dgCMatrix")),
                                    Rstruct = chol0))
    } else {
      
      if(is.null(loc_obs)){
        # Original sampler, based on Matrix package:
        
        # Cholesky of Quadratic term:
        chQht_Matrix = Matrix::chol(QHt_Matrix)
        
        # Sample the states:
        mu = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
        
      }else{
        if (D == 1) {
          diag1 = 1/obs_sigma_t2 + 1/evol_sigma_t2 + c(1/evol_sigma_t2[-1], 0)
          diag2 = -1/evol_sigma_t2[-1]
          rd = RcppZiggurat::zrnorm(T)
          mu = as.matrix(sample_mat_c(loc_obs$r, loc_obs$c, c(diag1, diag2, diag2), length(diag1), length(loc_obs$r), c(linht), rd, D))
        } else {
          diag1 = 1/obs_sigma_t2 + 1/evol_sigma_t2 + c(0, 4/evol_sigma_t2[-(1:2)], 0) + c(1/evol_sigma_t2[-(1:2)], 0, 0)
          diag2 = c(-2/evol_sigma_t2[3], -2*(1/evol_sigma_t2[-(1:2)] + c(1/evol_sigma_t2[-(1:3)],0)))
          diag3 = 1/evol_sigma_t2[-(1:2)]
          rd = RcppZiggurat::zrnorm(T)
          mu = as.matrix(sample_mat_c(loc_obs$r, loc_obs$c, c(diag1, diag2, diag2, diag3, diag3), length(diag1), length(loc_obs$r), c(linht), rd, D))
        }
      }
    }
  }
  
  # And return the states:
  mu
}
sampleEvol0 = function(mu0, evolParams0, commonSD = FALSE, A = 1){
  
  # Store length locally:
  p = length(mu0)
  
  # For numerical stability:
  mu02offset = any(mu0^2 < 10^-16)*max(10^-8, mad(mu0)/10^6)
  mu02 = mu0^2 + mu02offset
  
  if(commonSD){
    # (Common) standard deviations:
    evolParams0$sigma_w0 = rep(1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(mu02)/2 + evolParams0$px_sigma_w0[1])), p)
    
    # (Common) paramater expansion:
    evolParams0$px_sigma_w0 = rep(rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0[1]^2 + 1/A^2), p)
    
  } else {
    # (Distinct) standard deviations:
    evolParams0$sigma_w0 = 1/sqrt(rgamma(n = p, shape = 1/2 + 1/2, rate = mu02/2 + evolParams0$px_sigma_w0))
    
    # (Distinct) paramater expansion:
    #evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/A^2)
    evolParams0$px_sigma_w0 = rgamma(n = p, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_w0^2 + 1/evolParams0$sigma_00^2)
    
    # Global standard deviations:
    evolParams0$sigma_00 = 1/sqrt(rgamma(n = 1, shape = p/2 + 1/2, rate = sum(evolParams0$px_sigma_w0) + evolParams0$px_sigma_00))
    
    # (Global) parameter expansion:
    evolParams0$px_sigma_00 = rgamma(n = 1, shape = 1/2 + 1/2, rate = 1/evolParams0$sigma_00^2 + 1/A^2)
  }
  
  # And return the list:
  evolParams0
}
initEvolParams_HS <- function(omega, Td) {
  tau = 1
  tau_x = extraDistr::rinvgamma(1, 1, 1 + 1 / tau ^ 2)
  lambda_x = rep(1, Td)
  lambda = sqrt(extraDistr::rinvgamma(Td, 1, 1 / lambda_x + (omega/tau)^2/2))
  return(list(
    sigma_wt = pmax(tau * lambda, 1e-7),
    tau = tau,
    tau_x = tau_x,
    lambda = lambda,
    lambda_x = lambda_x
  ))
}
sampleEvolParams_HS <- function(omega, evolParams, Td) {
  tau = evolParams$tau
  tau_x = evolParams$tau_x
  lambda = evolParams$lambda
  lambda_x = evolParams$lambda_x
  
  tau = sqrt(extraDistr::rinvgamma(1, (Td + 1) / 2, 1 / tau_x + sum((omega/
                                                                       lambda) ^ 2) / 2))
  tau_x = extraDistr::rinvgamma(1, 1, 1 + 1 / tau ^ 2)
  lambda = sqrt(extraDistr::rinvgamma(Td, 1, 1 / lambda_x + (omega/tau)^2/2 ))
  lambda_x = extraDistr::rinvgamma(Td, 1, 1 + 1 / lambda ^ 2)
  
  evolParams$sigma_wt = tau * lambda
  evolParams$tau = tau
  evolParams$tau_x = tau_x
  evolParams$lambda = lambda
  evolParams$lambda_x = lambda_x
  return(evolParams)
}
initEvolParams_NIG <- function(omega,Td) {
  return(list(sigma_wt = tcrossprod(rep(1,Td), apply(omega, 2, function(x) sd(x, na.rm=TRUE)))))
}
sampleEvolParams_NIG <- function(omega,evolParams,Td){
  evolParams = list(sigma_wt = tcrossprod(rep(1,Td),
                                          apply(omega, 2,
                                                function(x) 1/sqrt(rgamma(n = 1, 
                                                                          shape = 0.01+Td/2, 
                                                                          rate = 0.01+sum(x^2)/2)))))
  return(evolParams)
}
initEvolParams_HS_sparse <- function(omega,
                                     Td,
                                     tau = 1 / (1000 * Td)) {
  #nmin1 = 0
  omega_norm = omega / tau
  x_lambda_t = rep(100, Td)
  lambda_2 = extraDistr::rinvgamma(Td, 1, 1 / x_lambda_t + omega_norm ^
                                     2 / 2)
  sigma_wt = pmax(sqrt(lambda_2) * tau, 1e-8)
  return(list(
    sigma_wt = sigma_wt,
    tau = tau,
    lambda_2 = lambda_2
  ))
}
sampleEvolParams_HS_sparse <- function(omega,
                                       Td,
                                       evolParams,
                                       tau = 1 / (1000 * Td)){
  tau = evolParams$tau
  lambda_2 = evolParams$lambda_2
  omega_norm = omega / tau
  x_lambda_t = extraDistr::rinvgamma(Td, 1, 1 + 1 / lambda_2)
  lambda_2 = extraDistr::rinvgamma(Td, 1, 1 / x_lambda_t + omega_norm ^
                                     2 / 2)
  sigma_wt = pmax(sqrt(lambda_2) * tau, 1e-8)
  return(list(
    sigma_wt = sigma_wt,
    tau = tau,
    lambda_2 = lambda_2
  ))
}
init_HSplus = function(omega,Td){
  xi = extraDistr::rinvgamma(1, 1/2, 1)
  tau_t2 = extraDistr::rinvgamma(1, 1/2, 1/xi)
  
  v = extraDistr::rinvgamma(Td, 1/2, 1/tau_t2)
  lambda_t2 = extraDistr::rinvgamma(Td, 1/2, 1/v)
  
  list(xi=xi, v=v, tau_t2 = tau_t2, lambda_t2 = lambda_t2, sigma_wt = sqrt(lambda_t2))
}
sample_HSplus = function(omega, evolParams,Td){
  #omega = as.matrix(omega)
  
  hsInput2 = omega^2
  
  evolParams$lambda_t2 = extraDistr::rinvgamma(Td, 1, 1/evolParams$v + hsInput2/2)
  evolParams$v = extraDistr::rinvgamma(Td, 1, 1/evolParams$lambda_t2 + 1/evolParams$tau_t2)
  
  evolParams$tau_t2 = extraDistr::rinvgamma(1, (Td+1)/2, 1/evolParams$xi + sum(1/evolParams$v))
  evolParams$xi = extraDistr::rinvgamma(1, 1, 1+1/evolParams$tau_t2)
  
  
  evolParams$sigma_wt = sqrt(evolParams$lambda_t2)
  
  evolParams
}

fit_DIRT = function(y, D = 1,evol_error = "HS",
                    nsave = 1000, nburn = 1000, nskip = 4,
                    verbose = TRUE){
  dims = dim(y)
  n_students = dims[1]
  n_questions = dims[2]
  N = dims[3]
  Td = N - D
  
  param_list = list()
  # intializing parameters
  for(i in 1:n_students){
    # transformed data
    kappa_i = colSums(y[i,,])-n_questions/2
    # precision parameter for the Conditional Gaussian likelihood
    w = rep(1,N)
    # variance of the first mu 
    evolParams0 = initEvol0(  kappa_i[1:D])
    # variance of the omega
    if(evol_error == "HS"){
      evolParams = initEvolParams_HS(  kappa_i[-c(1:D)],Td)    
    }else if(evol_error == "NIG"){
      evolParams = initEvolParams_NIG(as.matrix(kappa_i[-c(1:D)]),Td)    
    }else{
      evolParams = initEvolParams_HS_sparse(kappa_i[-c(1:D)],Td,tau = 1/Td)
    }
    
    # Student Proficiency
    theta = sampleBTF(kappa_i/w,
                      obs_sigma_t2 = 1/w,
                      evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                        evolParams$sigma_wt^2),
                      D = D)
    # difference of theta
    omega = diff(theta, differences = D)  
    theta0 = as.matrix(theta[1:D,])
    param_list[[i]] = list(kappa_i =   kappa_i)
    param_list[[i]]$w = w
    param_list[[i]]$evolParams0 = evolParams0
    param_list[[i]]$evolParams = evolParams
    param_list[[i]]$theta = theta
    param_list[[i]]$theta0 = theta0
    param_list[[i]]$omega = omega
  }
  #
  
  # Store the MCMC output in separate arrays (better computation times)
  mcmc_output = vector('list', 5); 
  names(mcmc_output) = c("theta","p","evol_sigma_t2","tau","lambda")
  post_theta = array(NA, c(nsave,n_students,N))
  post_p = array(NA,c(nsave,n_students,N))
  post_evol_sigma_t2 = array(NA, c(nsave,n_students, N))
  post_tau  = array(NA,c(nsave,n_students))
  post_lambda = array(NA, c(nsave,n_students,Td))
  
  # Total number of MCMC simulations:
  nstot = nburn+(nskip+1)*(nsave)
  skipcount = 0; isave = 0 # For counting
  
  # Run the MCMC:
  if(verbose){
    pb <- progress::progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Estimated time remaining: :eta]",
                                     total = nstot,
                                     complete = "=",   # Completion bar character
                                     incomplete = "-", # Incomplete bar character
                                     current = ">",    # Current bar character
                                     clear = FALSE,    # If TRUE, clears the bar when finish
                                     width = 100)      # Width of the progress bar
  }
  for(nsi in 1:nstot){
    if(verbose){
      if(nsi < 10){
        pb$tick()
      }
      else if(((nsi%%1000) == 0)){
        pb$tick(1000)
      }
    }
    
    # Sampling
    for(i in 1:n_students){
      kappa_i= param_list[[i]]$kappa_i
      evolParams0 = param_list[[i]]$evolParams0
      evolParams = param_list[[i]]$evolParams
      w = param_list[[i]]$w

      theta = sampleBTF(kappa_i/w,
                        obs_sigma_t2 = 1/w,
                        evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                          evolParams$sigma_wt^2),
                        D = D)
      w = pgdraw::pgdraw(n_questions,theta)
      theta0 = theta[1:D]
      omega = diff(theta,differences = D)
      evolParams0 = sampleEvol0(theta0,evolParams0)
      
      if(evol_error == "HS"){
        evolParams = sampleEvolParams_HS(omega,evolParams = evolParams,Td)
      }else if(evol_error == "NIG"){
        evolParams = sampleEvolParams_NIG(as.matrix(omega),evolParams = evolParams,Td)    
      }else{
        evolParams = sampleEvolParams_HS_sparse(omega,Td = Td,evolParams = evolParams,tau = 1/Td)    
      }
      
      param_list[[i]]$w = w
      param_list[[i]]$evolParams0 = evolParams0
      param_list[[i]]$evolParams = evolParams
      param_list[[i]]$theta = theta
      param_list[[i]]$theta0 = theta0
      param_list[[i]]$omega = omega
    }
    
    # Store the MCMC output:
    if(nsi > nburn){
      # Increment the skip counter:
      skipcount = skipcount + 1
      
      # Save the iteration:
      if(skipcount > nskip){
        # Increment the save index
        isave = isave + 1
        
        # Save the MCMC samples:
        for(i in 1:n_students){
          post_theta[isave,i,] = param_list[[i]]$theta
          post_p[isave,i,] = exp(param_list[[i]]$theta)/(1+exp(param_list[[i]]$theta))
          post_evol_sigma_t2[isave,i,] = c(param_list[[i]]$evolParams0$sigma_w0^2,
                                           param_list[[i]]$evolParams$sigma_wt^2)
          if(evol_error != "NIG"){
            post_lambda[isave,i,] = param_list[[i]]$evolParams$lambda
            post_tau[isave,i] = param_list[[i]]$evolParams$tau 
          }
        }
        
        
        
        # And reset the skip counter:
        skipcount = 0
      }
    }
  }
  
  mcmc_output$theta = post_theta
  mcmc_output$p = post_p
  mcmc_output$evol_sigma_t2 = post_evol_sigma_t2
  if(evol_error != "NIG"){
    mcmc_output$lambda = post_lambda
    mcmc_output$tau = post_tau
  }
  return (mcmc_output);
}



