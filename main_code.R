# measles
# measles escape probability
#############################################################################
#
# SCRIPT for homogeneous case: 
# q (escape probability) as same across all households
#
#############################################################################
# ---------------------------------------------------------------------------
# This code provides a function to draw samples from the joint posterior 
# distribution of the escape probability q and the frequency
# n111 of chain 1->1->1, based on the NYC measles 
# outbreak data. The algorithm is based on Gibbs sampling.
#
# INPUT mcmc.size  = the number of MCMC iterations
#       alpha,beta = parameters of the Beta prior distribution
#                    of the escape probability
#
# OUTPUT chainGibbs = the MCMC sample of q and n111 as a list
#                     of two entries, each a vector of length mcmc.size
# 
# Code v.1 on July 22, 2015. Updates on July 29-31, 2015
# Last update August 5, 2015
# ----------------------------------------------------------------------------
chainGibbs = function(mcmc.size, alpha, beta){
# Reserve space
q    = rep(0, mcmc.size)
p    = rep(0, mcmc.size)
n111 = rep(0, mcmc.size)

# Initialize the model unknowns
q[1]    = 0.5;                         # just an initial guess
p[1]    = 1 - q[1]
n111[1] = round(275*2*q[1]/(2*q[1]+1)) # the (rounded) expected value of n_111, 
                                       # given q = 0.5 

# The observations (see section 2 for details)
n1  = 34   # frequency of chain 1     (frequency is observed)
n11 = 25   # frequency of chain 1->1  (frequency is observed)
N3  = 275  # frequency of chains with outbreak size 3 (not observed)
           # (the total frequency of chains 1->2 and 1->1->1)

# Draw MCMC samples
for (i in 2:mcmc.size){

   # Draw an iterate of q from its full conditional distribution
   q[i] = rbeta(1,2*n1+2*n11+n111[i-1]+alpha,n11+2*N3+beta)
   p[i] = 1 - q[i]
  
   # Draw an iterate of n111 from its full conditional distribution
   n111[i] = rbinom(1,N3,2*q[i]/(2*q[i]+1))
   
   }

# The output: the MCMC samples
chainGibbs = list(q=q,n111=n111,p=p)

}

# mcmc-size = 5000, alpha = 1, beta = 1
test=chainGibbs(5000,1,1)

# Set up to have 2 plots in one figure
par(mfrow=c(1,2),oma=c(0,0,0,0))

# Assume a 'burn-in' period of 500 iterations, only plotting those greater than 500
# Plot iterates of q and n111 
hist(test$q[501:5000],main="",xlab="q")
hist(test$n111[501:5000],main="",xlab="n111")

# Summary of the results
summary(test$q[501:5000])
quantile(test$q[501:5000], c(0.025,0.975))
summary(test$n111[501:5000])
quantile(test$n111[501:5000], c(0.025, 0.975))

# ----------------------------------------------------------------------------
# checkmodel draws numerical samples from the joint posterior 
# predictive distribution of frequencies n1, n11, n111 and n12.
# Each individual sample is based on one realisation from the posterior
# distribution of the escape probability q. (The escape probability is
# common across all households.)
#
# The function outputs the mean of the predictive distribution, based on the 
# sample. It also plots the samples from the 
# (marginal) predictive distribution of frequencies n1 and n11. 
#
# INPUT mcmc.sample = the output (= MCMC samples) from chainGibbs
#       mcmc.burnin = the number of MCMC samples to be discarded
#
# OUTPUT expfreq = the posterior predictive expectations of the 
#                  four chain frequencies (n1,n11,n111,n12)
#
# Code v.1 on July 22, 2015. Updates on July 29-31, 2015
# Last update August 5, 2015
# ----------------------------------------------------------------------------
checkmodel = function(mcmc.sample,mcmc.burnin){

  # The number of MCMC samples in the input
  mcmc.size = length(mcmc.sample$q)
  
  # Discard the burn-in samples and retain the rest of the samples for 
  # the escape probability q
  q = mcmc.sample$q[(mcmc.burnin+1):mcmc.size]
  
  # The number of samples retained (= the size of the sub-sample)
  mcmc.subsamplesize = mcmc.size-mcmc.burnin
  
  # Reserve space for the predictive frequencies
  Npred = matrix(0,mcmc.subsamplesize,4)
  
  # Calculate the posterior predictive chain probabilities for
  # each retained MCMC sample. This produces vectors of length 
  # 'mcmc.subsamplesize': 
  P1   =  q^2
  P11  =  2*(1-q)*q^2   
  P111 =  2*(1-q)^2*q
  P12  = (1-q)^2
  
  # Sample from the posterior predictive distribution of chain frequencies
  # with the total  number of chains as 334.
  for (k in 1:mcmc.subsamplesize){
    Npred[k,] = rmultinom(1,334,c(P1[k],P11[k],P111[k],P12[k]))
  } 
  
  # Plot the posterior predictive distribution for frequencies n1 and n11
  par(mfrow=c(1,1))
  plot(0:60,0:60,type="n",xlab='n1',ylab='n11',cex.lab=2)           
  points(jitter(Npred[,1],2.5),jitter(Npred[,2],2.5),pch='.',cex=3) 
  
  # The actually observed value
  points(34,25,col='red')
  abline(v=34,col='red',lwd=0.25)
  abline(h=25,col='red',lwd=0.25)
  
  # The posterior predictive expected frequnecies of the four chain types
  expfreq = round(334*c(mean(P1),mean(P11),mean(P111),mean(P12)))
  
  # Return the posterior predictive expected frequencies
  print("The posterior predictive expected frequencies of the four chain types")
  expfreq
  
  return(expfreq)
  
}

# Running the functions: draw a posterior sample and check the model fit.
mcmc.sample = chainGibbs(5000,1,1)

checkmodel(mcmc.sample,mcmc.burnin=500)

#                      END of homogeneous case script  
############################################################################

############################################################################
# 
# SCRIPT heterogeneous case: hierarchical model 
# q (escape probability) varies across households
#
############################################################################

## -------------------------------------------------------------------------
## This code provides routines to draw samples from the joint posterior 
# distribution of parameters 'tildeq' and 'z' in a hierarchical chain 
# binomial model of the measles outbreak data. The model unknowns include
# the household-specific escape probabilities. 
#
# (n1^j,n11^j,n111^j,n12^j)| q_j ~ Multinom(q_j^2,2q_j^2p_j,2q_jp_j^2,p_j^2)
#  q_j| tildeq, z                ~ Beta(tildeq/z,(1-tildeq)/z)
#  tildeq                        ~ Uniform(0,1)
#  z                             ~ Gamma(1.5,1.5)
#
## Code v.1 on July 29, 2015. Updates on August 1-4, 2015
## Last update August 5, 2015
## -------------------------------------------------------------------------
propose.par = function(cur.par,prop.unif,delta.w){
  # This function samples a new (candidate) value  
  # around the current value of the parameter iterate
  return(cur.par + (0.5-prop.unif)*delta.w)
}

log.likelihood = function(my.tildeq,my.z,q.vec){
  # This function calculates the log-likelihood function of the two parameters,
  # based on the (current iterates) of the household-specific 
  # escape probabilities
  
  # INPUT my.tildeq = parameter tildeq
  #       my.z.     = parameter z
  #       q.vec     = a vector of household-specific escape probabilities
  
  return(sum(log(dbeta(q.vec,my.tildeq/my.z,(1-my.tildeq)/my.z))))
}

chain_hierarchical = function(mcmc.size){
  # This is the main sampling function. The data, the parameters of the 
  # prior distributions and the 'tuning parameters' (i.e the proposal window widths)
  # are hard-code. 
  # In addition to the two model parameters ('tildeq' and z), the 
  # algorithm draws samples of the 334 escape probabilities and 
  # the unknown realisations of chains in the 275 households with 
  # the final number infected as 3.
  # INPUT mcmc.size = the number of MCMC iterations to be sampled
  #
  # OUTPUT chain_hierarchical = the MCMC sample of parameters 'tildeq' and 'z'
  #                             of the Beta distribution for household specific
  #                             escape probabilities
  
  # Parameters of the Gamma(nu1,nu2) prior for parameter z
  nu1 = 1.5
  nu2 = 1.5 
  # Tuning parameters: the widths of the proposal windows
  delta.tildeq = 0.08 # for parameter q
  delta.z      = 0.3  # for parameter z
  
  # Step (1): Reserve space for the MCMC samples
  tildeq  = rep(0,mcmc.size)
  z       = rep(0,mcmc.size)
  q       = matrix(0,nrow=mcmc.size,ncol=334)
  n111    = matrix(0,nrow=mcmc.size,ncol=275)
  
  # Step (2): Initialize variables
  cur.tildeq = 0.5
  cur.z      = 1.0
  tildeq[1]  = cur.tildeq
  z[1]       = cur.z
  q[1,]      = rep(0.5,334)
  n111[1,]   = rbinom(1,275,2*0.5/(2*0.5+1))
  
  # The observations
  n1  = 34   # frequency of chain 1
  n11 = 25   # frequency of chain 1->1
  N3  = 275  # frequency of chains of size 3 (i.e., 1->2 or 1->1->1)
  
  # Draw MCMC samples
  for (i in 2:mcmc.size){
    
    # Step (3): Update the household-specific escape probabilities
    #           via Gibbs procedure
    
    alpha =  cur.tildeq/cur.z     # tildeq[i-1]/z[i-1]
    beta  =  (1-cur.tildeq)/cur.z #(1-tildeq[i-1])/z[i-1]
    
    # households with chain 1
    for (j in 1:n1){
      q[i,j] = rbeta(1,2+ alpha,beta)
    }
    
    # households with chain 1->1
    for (j in 1:n11){
      q[i,j+n1] = rbeta(1,2+alpha,1+beta)
    }
    
    # households with chain 1->1->1 or 1->2
    for (j in 1:275){
      q[i,j+n1+n11] = rbeta(1,n111[i-1,j]+alpha,2+beta)
    }
    
    # Step (4) Update the household-specific frequencies of chain 1->1->1
    #          via Gibbs procedure

    for (j in 1:N3){ 
      
      # Draw a new iterate of n111 from its full conditional distribution
      n111[i,j] = rbinom(1,1,2*q[i,j+n1+n11]/(2*q[i,j+n1+n11]+1))
    }
    
    # Step (5) Update parameter tildeq via Metropolis-Hastings 
    # Propose a new value for parameter tildeq
    new.tildeq = propose.par(cur.tildeq,runif(1),delta.tildeq)
    
    # If the proposed value if within the range of the parameter 
    if ((new.tildeq >0) & (new.tildeq < 1)){
      
      # The log acceptance ratio = the log likelihood ratio   
      #(a uniform prior is assumed for tildeq and a symmetric proposal has been made)
      log.acc.tildeq = log.likelihood(new.tildeq,cur.z,q[i,]) - log.likelihood(cur.tildeq,cur.z,q[i,]) 
      
      # Metropolis step for tildeq: determine whether the proposal is accepted
      if (log(runif(1)) < log.acc.tildeq){
        cur.tildeq = new.tildeq
      }
    }
    
    # Step (6): Update parameter z via Metropolis-Hastings
    # Propose a new value for parameter z
    new.z  = propose.par(cur.z,runif(1),delta.z)
    
    # If the proposed value is withing the range of the parameter
    if (new.z > 0) { 
      
      # The log likelihood ratio
      log.like.z = log.likelihood(cur.tildeq,new.z,q[i,]) - log.likelihood(cur.tildeq,cur.z,q[i,]) 
      
      # The log acceptance ratio = the log-posterior ratio (since a symmetric proposals are being made)
      log.acc.z = log.like.z + log(dgamma(new.z,nu1,nu2)) - log(dgamma(cur.z,nu1,nu2))
      
      # Metropolis step for z: determine whether the proposal is accepted
      if (log(runif(1))< log.acc.z){
        cur.z = new.z
      }
    }   
    # Save the current iterates of parameters tildeq and z
    tildeq[i] = cur.tildeq
    z[i]      = cur.z 
  }
  # The output: the MCMC samples
  return(list(tildeq=tildeq,z=z))
}

# Running the script to draw a posterior sample of size 2000 
mcmc.sample = chain_hierarchical(mcmc.size=2000)

# Plot the posterior of tilde q
hist(mcmc.sample$tildeq[500:2000],xlab='tilde q',xlim=c(0.1,0.35))

# -----------------------------------------------------------------
# Script for checking the model
# -----------------------------
check_hierarchical = function(mcmc.sample,mcmc.burnin){
  #
  # This function draws samples from the posterior predictive 
  # distribution of chain frequencies (n1, n11, n111, n12)
  # in the hierarchical model.
  #
  # The function plots the samples from the (marginal) 
  # distribution of frequencies (n1,n11) for comparison with the
  # actually observed value.
  #
  # The fuction output is the posterior expectations of the four
  # chain frequencies in a sample of 334 households
  #
  # INPUT mcmc.sample = the output from program chain_hierarchical.R
  #       mcmc.burnin = the number of samples to be discarded from the MCMC sample
  #
  # OUTPUT check_hierarchical = the posterior predictive expectations of the 
  #                             four frequencies
  
  # The number of MCMC samples
  mcmc.size = length(mcmc.sample$tildeq)
  
  # Discard the burn-in samples
  tildeq = mcmc.sample$tildeq[(mcmc.burnin+1):mcmc.size]
  z      = mcmc.sample$z[(mcmc.burnin+1):mcmc.size]
  
  # The number of retained samples
  mcmc.subsamplesize = mcmc.size-mcmc.burnin
  
  # Reserve space for the predictive samples of the chain frequencies
  Npred = matrix(0,mcmc.subsamplesize,4)
  
  # Reserve space for one chain realisation in each of the 334 households  
  L    = matrix(0,334,4) 
  
  # Iterate over the MCMC samples
  for (i in 1:mcmc.subsamplesize){
    
    # Sample household-specific escape probabilities and
    # the chains in 334 households
    for (j in 1:334){
      
      q     = rbeta(1,tildeq[i]/z[i],(1-tildeq[i])/z[i])
      P1    = q*q
      P11   = 2*(1-q)*q^2
      P111  = 2*q*(1-q)^2
      P12   = (1-q)^2
    
      # A realisation of the chain
      L[j,] = rmultinom(1,1,c(P1,P11,P111,P12))
    }
    
    # Store the i'th MCMC sample of the chain frequencies
    Npred[i,] = apply(L,2,sum)
  } 
  
  # Plot the sample points from the posterior predictive distribution 
  # of frequencies n1 and n11
  plot(0:60,0:60,type="n",xlab='n1',ylab='n11',cex.lab=2)
  points(jitter(Npred[,1],2.5),jitter(Npred[,2],2.5),pch='.',cex=3)
  
  # The actually observed value (n1,n11)
  points(34,25,col='red')
  abline(v=34,col='red',lwd=0.25)
  abline(h=25,col='red',lwd=0.25)
  
  # Return the posterior predictive expected frequencies
  check_hierarchical = return(round(apply(Npred,2,mean)))
  
}

# Checking the fit of the hierarchical model
check_hierarchical(mcmc.sample,mcmc.burnin=500)

#                     END of hierarchical model script
##############################################################################


