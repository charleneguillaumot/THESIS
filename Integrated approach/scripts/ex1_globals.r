# ex1_globals.r
# 
#   Copyright 2014 Matthew V Talluto, Isabelle Boulangeat, Dominique Gravel
# 
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 3 of the License, or (at
#   your option) any later version.
#   
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#   General Public License for more details.
# 
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
# 
# 
# settings and functions used in multiple scripts for example 1


#library(pacman)
library(rjags)
library(pbapply)

settings <- list(
	chains=2, 
	tuning = 1000, 
	burnin = 4000, 
	samples=2000, 
	thin=10, 
	seed=NA,
	cv_samples=100,
	diagnostics=TRUE,
	quantiles=c(0.05, 0.95))

## information for model 1 and the metamodel -- where to map temp and precipitation
tempLevels <- seq (6.62, 7.36, 0.1)
fLevels <- seq(0,1, 0.1)	
depthLevels <- seq(-195,0, 0.1)	
# levels of temperature, f and depth for predictions
predictionMap <- expand.grid(tempLevels, fLevels,depthLevels) # data frame for every possible combination of T and P
names(predictionMap) <- c("temp", "f", "depth")

#precipRegimes <- list(current=c(0.1, 1.1), future=c(-1.1, 0.3)) # current and future range for precipitation
#precipPrediction <- seq(min(unlist(precipRegimes)), max(unlist(precipRegimes)), length.out=50) # prediction levels for precipitation only predictions



###############
# helper functions
###############
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

logit_inv <- function(x) {
  inf <- which(is.infinite(x) & x > 0)
  ninf <- which(is.infinite(x) & x < 0)
  p <- exp(x) / (1 + exp(x))
  p[inf] <- 1
  p[ninf] <- 0
  return(p)
}

prediction_summary <- function(p, quantiles=c(0.05, 0.95)) {
  M <- sapply(1:ncol(p), function(i) mean(p[,i],na.rm=T))
  SD <- sapply(1:ncol(p), function(i) sd(p[,i],na.rm=T))
  Q1 <- sapply(1:ncol(p), function(i) quantile(p[,i], min(quantiles),na.rm=T))
  Q2 <- sapply(1:ncol(p), function(i) quantile(p[,i], max(quantiles),na.rm=T))
  
  return(data.frame(mean=M, SD=SD, lower=Q1, upper=Q2))
}

predict_psi <- function(posterior, var1, var2 = NULL, summary=TRUE, logit=TRUE, ...) {
  M <- mcmc_unlist(posterior)
  predictions <- t(sapply(1:nrow(M), function(i) calculate_psi(M[i,], var1, var2, logit)))
  if(summary) predictions <- prediction_summary(predictions, ...)
  
  return(predictions)
}

calculate_psi <- function(pars, x1, x2 = NULL, logit=TRUE) {
  yu <- pars[1] + pars[2]*x1 + pars[3]*x1^2
  if(!is.null(x2))
    yu <- yu + pars[4]*x2 + pars[5] * x2^2
  if(logit) yu <- logit_inv(yu)
  return(yu)
}

mcmc_unlist <- function(M) {
  r <- M[[1]]
  if(length(M) > 1) {
    for(i in 2:length(M))
      r <- rbind(r, M[[i]])
  }
  return(r)
}



formula_poly_apply <- function(x,b){#apply formula_poly for all rows in x
  apply(x,1,function(y)formula_poly(y,b)) 
}

predictions_poly <- function(x,b){#to call
  return(t(pbapply(b,1,function(y)formula_poly_apply(x,y))))
}


formula_poly <- function(x,b){
  
         boot::inv.logit(b[1] + b[2]*x[1] + 
                           b[3] * x[1]^2+ b[4] * x[2] + b[5] * x[2]^2)

}



formula <- function(x,b){
  #if (aicc_all[1]==aicc_poly){
    boot::inv.logit(b[1] + b[2]*x[1] + 
                      b[3] * x[2]+ b[4] * x[3] + b[5] * x[2]^2) + b[6] * x[3]^2}#}

formula_apply <- function(x,b){
  apply(x,1,function(y)formula(y,b)) 
}

predictions <- function(x,b){
  return(t(pbapply(b,1,function(y)formula_apply(x,y))))
}




# 
predictions_poly <- function(x,b){
  prob_data <- data.frame(matrix(NA, nrow = nrow(b), ncol = nrow(x)))
  for (i in 1:nrow(b)){

    prob_data[i,] <- apply(x,1,formula_poly,b[i,])
    print(i)

  }
  return(prob_data)
}


# predictions_poly_asp_pad <- function(x,b){
#   prob_data <- data.frame(matrix(NA, nrow = nrow(b), ncol = nrow(x)))
#   for (i in 1:nrow(b)){
#     
#     prob_data[i,] <- apply(x,1,formula_poly_asp_pad,b[i,])
#     print(i)
#     
#   }
#   return(prob_data)
# }

ploting <- function(m,title,xlabel){
  prob_sstmax_plot <- ggplot(m, aes(x=m[,1], y=m[,2])) + 
#     geom_errorbar(aes(ymin=(m[,2])-(m[,3]), ymax=(m[,2])+(m[,3])),
#                   width=.3, color = "black") +
    geom_line(size=1, color = "blue") +
    geom_point(color = "black")+
    geom_ribbon(aes(ymin = m[,4], ymax = m[,5], fill="Uncertainty"),alpha = .25)
  prob_sstmax_plot <- prob_sstmax_plot + labs(x=(xlabel), y="Probability")
  prob_sstmax_plot <- prob_sstmax_plot + ggtitle(title)
  prob_sstmax_plot
}


## generate the model files
# cat("
# data {
#     N <- length(presence)
# }
# model {
# 
# for (i in 1:N) {
#     presence[i] ~ dbin(pr[i], 1)
#     logit(pr[i]) <- b0 + b1 * temp[i] + b2 * temp[i]^2 + b3 * precip[i] + b4 * precip[i]^2
# }
# 
# # priors
# b0 ~ dnorm(bPrior[1,1], bPrior[1,2])
# b1 ~ dnorm(bPrior[2,1], bPrior[2,2])
# b2 ~ dnorm(bPrior[3,1], bPrior[3,2])
# b3 ~ dnorm(bPrior[4,1], bPrior[4,2])
# b4 ~ dnorm(bPrior[5,1], bPrior[5,2])
# }
# "
# , file="ex1_metamodel.jags")
# 
# 
# cat("
# data {
# 	N <- length(psi)
# }
# model {
# for (i in 1:N) {
# 	psi[i] ~ dbeta(p[i], q[i])
# 	p[i] <- mu[i] * phi
# 	q[i] <- (1-mu[i]) * phi
# 	
# 	logit(mu[i]) <- a0 + a1*precip[i] + a2 * precip[i]^2
# }
# 
# # priors
# a0 ~ dnorm(0, 0.001)
# a1 ~ dnorm(0, 0.001)
# a2 ~ dnorm(0, 0.001)
# 
# ## uninformative prior for phi from Gelman 2006
# phi <- U^2
# U ~ dunif(0,50)
# }"
# , file="ex1_model2.jags")
# 
# cat("
# data {
# 	N <- length(psi)
# }
# model {
# for (i in 1:N) {
# 	psi[i] ~ dbeta(p[i], q[i])
# 	p[i] <- mu[i] * phi
# 	q[i] <- (1-mu[i]) * phi
# 	
# 	logit(mu[i]) <- a0 + a1*temperature[i] + a2 * temperature[i]^2
# }
# 
# # priors
#   a0 ~ dnorm(bPrior[1,1], bPrior[1,2])
#   a1 ~ dnorm(bPrior[2,1], bPrior[2,2])
#   a2 ~ dnorm(bPrior[3,1], bPrior[3,2])
# 
# ## uninformative prior for phi from Gelman 2006
# phi <- U^2
# U ~ dunif(0,50)
# }"
# , file = "phy_model2.jags")

### Globals for dealing with figures
## global variables
width <- 6.5	# dims of graphics window

# ## some settings for plots
title.line <- 0.7
label.line <- 1.5
axis.line <- 0.4

cex.title <- 0.75
cex.axis <- 0.7
lwd.axis <- 0.7
cex.lab <- 0.8
margins <- c(3,2.5,2.5,0.5)

bty='n'
mgp=c(label.line,axis.line, 0)
tcl=-0.3


make_axes <- function(xlab, ylab, main, ...) {
	axis(side=1, cex.axis=cex.axis, lwd=lwd.axis, lwd.ticks=lwd.axis, ...)
	axis(side=2, cex.axis=cex.axis, lwd=lwd.axis, lwd.ticks=lwd.axis, ...)
	title(xlab=xlab, ylab=ylab, cex.lab = cex.lab)
	mtext(main, side=3, line=title.line, cex=cex.title, adj=0)
}

