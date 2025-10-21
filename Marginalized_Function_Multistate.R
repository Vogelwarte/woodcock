#########################################
####By Jaume A. Badia-Boher##############
####Created December 2024################
#########################################


#Marginalized likelihood functions for multistate models
#Modified to allow for different states or departure
#The model should also work for cases where the state of departure is always the same. Recommended as a general model.
#The model does not treat uncertainty in state and event assignment at first encounter, hence it is exclusively for multistate models.

dmslik3 <- nimbleFunction(
  run = function(x = double(1), #Every function must start with argument z. Matrix y_sum which each unique capture history
                 
                 sumf = double(0), #First encounter for each unique capture history
                 nint = double(0), #Number of capture intervals, is nyears - 1
                 nstates = double(0), #Number of states
                 FR = double(0), #Vector with the number individuals that share every unique capture history
                 ps = double(3), #State transition matrix
                 po = double(3), #Observation transition matrix
                 log = integer(0, default = 0) #Compulsory argument, 0 for likelihood, 1 for loglikelihood
  ){ 
    
    returnType(double(0))
    zeta <- matrix(value = NA, nrow = nint+1, ncol = nstates)
    
    #We can't condition on first capture, as some individuals start at a 1, some at a 2.
    
    dvec <- numeric(length = nstates, value = 0)
    dvec[x[sumf]] <- 1
    
    zeta[sumf, 1:nstates] <- dvec
    
    #Calculate likelihood after first capture and loop by nstates as well
    
    for(t in sumf:nint){ 
      for(s in 1:nstates){
        
        zeta[t+1, s]  <- inprod(zeta[t,1:nstates], ps[1:nstates,s,t]) * po[s,x[t+1],t]
        
      }
    }
    
    
    lik <- sum(zeta[nint+1,1:nstates])
    logProb <- FR*log(lik)
    
    if(log) return(logProb)
    else return(exp(logProb))
    
  })

rmslik3 <- nimbleFunction(
  run = function(n = integer(),
                 sumf = double(0),
                 nint = double(0),
                 nstates = double(0),
                 FR = double(0), #Parameter not used, but needed as it has to match the params in the distribution function
                 ps = double(3),
                 po = double(3)
  ){
    
    returnType(double(1))
    z <- numeric(nint+1)
    
    #Individual at t=f can either be a 1 or a 2
    z[sumf] <- rcat(n = 1, prob = c(0.5,0.5))
    
    y <- z
    for(t in (sumf+1):(nint+1)){
      
      # state at t given state at t-1
      z[t] <- rcat(n = 1, prob = ps[z[t-1],1:nstates,t-1]) 
      # observation at t given state at t
      y[t] <- rcat(n = 1, prob = po[z[t],1:nstates, t-1]) 
      
    }
    
    return(y)
    
  }
)