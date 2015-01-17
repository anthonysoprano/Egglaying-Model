########################################################################################################
#      This script solves the egg-laying dynamics in C.elegans and fits the parameters so             ##
#      that the residual sum of squares is minimized w.r.t the experimental data                      ##
#                                                                                                     ##
#      QTL mapping is then performed using these parameters and resiudals as the phenotypes           ##
#      on the C.elegans genome                                                                        ##
########################################################################################################

# Load the required libraries 
library(deSolve)        
library(minpack.lm)
library(qtl)

#####################################       Egg-laying models        ###################################
## Egg-laying models are called by the ode functions. Given a list of times(t), parameters(parVar),   ##
## and initial states(stateVar), these functions will return a list of State variables at each of the ##
## times.                                                                                             ##
########################################################################################################

# Model 1

egglayingModel1 = function(t,stateVar,parVar) 
        ## dE = ko*O*S
        ## dS = -dE
        ## dO = kg - ks*O - dE
        ## parVar = [ko,kg,ks]
{       
        with(as.list(c(stateVar,parVar)),
{
        dE = ko*O*S
        dS = -dE
        dO = kg - ks*O - dE
        return(list(dS,dO,dE))
})
}

# Model 2 

egglayingModel2 = function(t,stateVar,parVar)
        ## dE = min(kf*O,ks*S)
        ## dS = -dE
        ## dO = ko - kc*O - dE
        ## parVar = [kf,ks,ko,kc]
{
        with(as.list(c(stateVar,parVar)),
{
        dE = min(kf*O,ks*S)
        dS = -dE
        dO = ko - kc*O - dE
        return(list(dS,dO,dE))
})
}

# Model 3 

egglayingModel3 = function(t,stateVar,parVar) 
        ## dE = kf*O*S
        ## dS = -dE
        ## dO = ko - dE
        ## parVar = [kg,ks,ko]
{
        with(as.list(stateVar,parVar),
{
        dE = kf*O*S
        dS = -dE
        dO = ko - dE
        return(list(dS,dO,dE))
})
}

# Model 4

egglayingModel4 = function(t,stateVar,parVar)
        ## dE = min(kf*O,ks*S)
        ## dS = -dE
        ## dO = ko,i - dE
        ## parVar = [kg,kf,ks]
{
        with(as.list(c(stateVar,parVar)),
{
        dE = min(kf*O,ks*S)
        dS = -dE
        dO = ko - dE
        return(list(dS,dO,dE))
})
}

######################################### The Fit Function ################################################