##Implementation of logistic regression with Polya-Gamma Data Augmentation
set.seed(150)
theta = -1
T = 200
theta_real = rep(c(0,2),each = 100) + theta
p = exp(theta_real)/(1+exp(theta_real))
y = rbinom(T,1,p)
# sampl size + meta data
burnin = 5000
samplesize = 5000
total = samplesize + burnin
theta_dat = array(NA,c(samplesize,T))

# initiaializing
D = 1
w = rep(1,T)
evolParams0 = dsp::initEvol0(y[1:D]-1/2)
evolParams = dsp::initEvolParams(y[-c(1:D)]-1/2,evol_error = "NIG")
for(i in 1:total){
  linht = y - 1/2
  QHt_Matrix = dsp::build_Q(obs_sigma_t2 = 1/w,evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                                                 evolParams$sigma_wt^2),
                            D = D)
  chQht_Matrix = Matrix::chol(QHt_Matrix)
  # Sample the states:
  theta = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
  w = pgdraw::pgdraw(1,theta)
  evolParams0 = dsp::sampleEvol0(theta[1:D],evolParams0)
  evolParams = dsp::sampleEvolParams(diff(theta,differences = D),evolParams = evolParams,evol_error = "NIG")
  if(i > burnin) theta_dat[i-burnin,] = theta
}

burnin = 20000
samplesize = 5000
total = samplesize + burnin
theta_hangsugraduation = array(NA,c(samplesize,T))
D = 1
w = rep(1,T)
evolParams0 = dsp::initEvol0(y[1:D]-1/2)
evolParams = dsp::initEvolParams(y[-c(1:D)]-1/2,evol_error = "HS")
for(i in 1:total){
  linht = y - 1/2
  QHt_Matrix = dsp::build_Q(obs_sigma_t2 = 1/w,evol_sigma_t2 = c(evolParams0$sigma_w0^2,
                                                                 evolParams$sigma_wt^2),
                            D = D)
  
  chQht_Matrix = Matrix::chol(QHt_Matrix)
  # Sample the states:
  theta = as.matrix(Matrix::solve(chQht_Matrix,Matrix::solve(Matrix::t(chQht_Matrix), linht) + rnorm(T)))
  w = pgdraw::pgdraw(1,theta)
  evolParams0 = dsp::sampleEvol0(theta[1:D],evolParams0)
  evolParams = dsp::sampleEvolParams(diff(theta,differences = D),evolParams = evolParams,evol_error = "HS")
  if(i > burnin) theta_hangsugraduation[i-burnin,] = theta
}

plot(theta_real,ylim = c(-4,4))
lines(colMeans(theta_hangsugraduation),col = "blue",lwd = 4)
lines(apply(theta_hangsugraduation,2,quantile,0.025),lwd =2,lty = 2,col = "blue")
lines(apply(theta_hangsugraduation,2,quantile,0.975),lwd =2,lty = 2,col = "blue")
lines(colMeans(theta_dat),col = "red", lwd = 4)
lines(apply(theta_dat,2,quantile,0.025),lwd =2,lty = 2,col = "red")
lines(apply(theta_dat,2,quantile,0.975),lwd =2,lty = 2,col = "red")

