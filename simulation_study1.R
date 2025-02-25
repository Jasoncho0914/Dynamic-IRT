#Code for Dynamic IRT with Time-Varying Covariate. 
# 10 students?
N = 30
theta1 = seq(from =-2,to = 2,length.out = N)
theta2 = rep(c(-1.5, 0, -1, 2, -1, 3),each = N/6)
theta3 = exp(theta1)/2
# item difficulty
nQ = 50
#item_difficulty = rnorm(30*nQ,mean = 0,sd = 0.5)
item_difficulty = rep(0,nQ * N)

data = array(NA,c(3,nQ,N))
for(t in 1:30){
  q_idx = (nQ*(t-1)+1):(nQ*t)
  dj = item_difficulty[q_idx]
  p1 = exp(theta1[t] - dj)/(1+exp(theta1[t] - dj))
  p2 = exp(theta2[t] - dj)/(1+exp(theta2[t] - dj))
  p3 = exp(theta3[t] - dj)/(1+exp(theta3[t] - dj))
  data[,,t] = rbind(rbinom(length(p1),1,p1),
                    rbinom(length(p2),1,p2),
                    rbinom(length(p3),1,p3))
}
#################################################################################
simulated_result_GL = fit_DIRT(data,D = 1,nsave = 5000,nburn = 5000,evol_error = "HS")
simulated_result_GL2 = fit_DIRT(data,D = 1,nsave = 5000,nburn = 5000,evol_error = "HSP") # more sparse than HS
simulated_result_NIG = fit_DIRT(data,D = 1,nsave = 5000,nburn = 5000,evol_error = "NIG")

plot(theta2,ylim = c(-3,4))
lines(colMeans(simulated_result_GL2$theta[,2,]),type = "l",col = "purple")
lines(apply(simulated_result_GL2$theta[,2,],2,quantile,0.025),type = "l",col = "purple",lty = 2)
lines(apply(simulated_result_GL2$theta[,2,],2,quantile,0.975),type = "l",col = "purple",lty = 2)


lines(colMeans(simulated_result_GL$theta[,2,]),type = "l",col = "blue")
lines(apply(simulated_result_GL$theta[,2,],2,quantile,0.025),type = "l",col = "blue",lty = 2)
lines(apply(simulated_result_GL$theta[,2,],2,quantile,0.975),type = "l",col = "blue",lty = 2)
lines(colMeans(simulated_result_NIG$theta[,2,]),type = "l",col = "red")
lines(apply(simulated_result_NIG$theta[,2,],2,quantile,0.025),type = "l",lty = 2,col = "red")
lines(apply(simulated_result_NIG$theta[,2,],2,quantile,0.975),type = "l",lty = 2,col = "red")

plot(theta1,ylim = c(-3,3))
lines(colMeans(simulated_result_GL$theta[,1,]),type = "l",col = "blue")
lines(apply(simulated_result_GL$theta[,1,],2,quantile,0.025),type = "l",col = "blue",lty = 2)
lines(apply(simulated_result_GL$theta[,1,],2,quantile,0.975),type = "l",col = "blue",lty = 2)
lines(colMeans(simulated_result_NIG$theta[,1,]),type = "l",col = "red")
lines(apply(simulated_result_NIG$theta[,1,],2,quantile,0.025),type = "l",lty = 2,col = "red")
lines(apply(simulated_result_NIG$theta[,1,],2,quantile,0.975),type = "l",lty = 2,col = "red")

plot(theta3,ylim = c(-3,4))
lines(colMeans(simulated_result_GL$theta[,3,]),type = "l",col = "blue")
lines(apply(simulated_result_GL$theta[,3,],2,quantile,0.025),type = "l",col = "blue",lty = 2)
lines(apply(simulated_result_GL$theta[,3,],2,quantile,0.975),type = "l",col = "blue",lty = 2)
lines(colMeans(simulated_result_NIG$theta[,3,]),type = "l",col = "red")
lines(apply(simulated_result_NIG$theta[,3,],2,quantile,0.025),type = "l",lty = 2,col = "red")
lines(apply(simulated_result_NIG$theta[,3,],2,quantile,0.975),type = "l",lty = 2,col = "red")
