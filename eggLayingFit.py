## Loading the required libraries
from lmfit import minimize, Parameters
from scipy.integrate import odeint
import matplotlib.pyplot as plt
import numpy as np  

def egglayingModel1(stateVar,t,params):
	## dE = ko*O*S
    ## dS = -dE
    ## dO = kg - ks*O - dE
	S = stateVar[0]
	O = stateVar[1]
	E = stateVar[2]
	kg = params[0]
	ks = params[1]
	ko = params[2]
	
	dE = ko*O*S
	dS = -dE
	dO = kg - ks*O -dE
	
	return [dS,dO,dE]


def egglayingModel2(t,state,params):
	
	## dE = min(kf*O,ks*S)
    ## dS = -dE
    ## dO = ko - kc*O - dE
    S = state[0]
    O = state[1]
    E = state[2]
    ko = params[0]
    kc = params[1]
    kf = params[2]
    ks = params[3]
    
    dE = min(kf*O,ks*S)
    dS = -dE
    dO = ko - kc*O - dE
    
    return [dS,dO,dE]

def egglayingModel3(t,state,params):
	## dE = kf*O*S
    ## dS = -dE
    ## dO = ko - dE
	S = state[0]
	O = state[1]
	E = state[3]
	ko = params[0]
	kf = params[1]
	
	dE = kf*O*S
	dS = -dE
	dO = ko - dE
	
	return [dS,dO,dE]	
   

def egglayingModel4(state,t,params):
	## dE = min(kf*O,ks*S)
	## dS = -dE
	## dO = ko - dE
	S = state[0]
	O = state[1]
	E = state[3]
	ko = params[0]
	kf = params[1]
	ks = params[2]
	
	dE = min(kf*O,ks*S)
	dS = -dE
	dO = ko - dE
	
	return [dS,dO,dE]
	
#def residual(params, state):

############################ Main Function #############################

## Parameter and State Initialization

## States
S0 = 300 # Initial No. of Sperms
O0 = 0	 # Initial No. of Oocytes
E0 = 0   # Initial No. of Eggs
state = [S0,O0,E0]

## Parameters Model1  
kg = 12
ks = 0
ko = 0.01
par1 = [kg,ks,ko]

## Parameters Model2
ko = 12
kc = 0
kf = 0.1
ks = 0.1
par2 = [ko,kc,kf,ks]

## Parameters Model3
ko = 12
kf = 0.01
par3 = [ko,kf]

## Parameters Model4
ko = 12
kf = 0.1
ks = 0.1
par4 = [ko,kf,ks]

## Time Gird
t = np.arange(0,96,1)
## Solve the ODEs
out = odeint(egglayingModel1, state, t,args = (par1,))
S = out[:,0]
O = out[:,1]
E = out[:,2]

## Plot the results
plt.figure()
plt.plot(t, S, label='Sperms')
plt.plot(t, O, label='Oocytes')
plt.plot(t, E, label='Eggs')
plt.title('Egglaying Dynamics-C.elegans')
plt.legend(loc=0)
plt.show()

