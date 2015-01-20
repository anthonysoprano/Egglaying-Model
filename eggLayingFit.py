## Loading the required libraries
from lmfit import minimize, Parameters
from scipy.integrate import odeint
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np  
import csv
import argparse

########################################################################
##						Command Line Options						  ## 
########################################################################
parser = argparse.ArgumentParser(description=
'This is a script to simulate and solve egglaying models in C.elegans')

parser.add_argument('-i','--input', help='Input file name for Experimental Data'
					,required=True)
parser.add_argument('-m','--model',help='Model No.(integer)(1/2/3/4)'
					, required=True)

args = parser.parse_args()
print "Experimental Data File: %s" % args.input 
print "Model No. Used: %s" % args.model 

########################################################################
##          Egglaying Models, take in globals containing state and    ##
##			the parameters to the model equations					  ##										
########################################################################

## Model1
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

## Model2
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

## Model3
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
   
## Model4
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
	
#def residual(params, model):
#	out = odeint(model, state, t,args = (params,))
#	modEggs = (out[round(e_times/delta_t), paste("E",i,sep="")]-out[round(e_times/delta_t)-1, paste("E",i,sep="")])/delta_t


########################################################################
##							 MAIN FUNCTION                            ##
########################################################################

## Experimental File Input
f = open('/home/raghav/git/Egglaying/DataSets/Egglaying_Model.csv','rU')
csv_f = csv.reader(f)


## Parameter and State Initialization

## States
S0 = 300 # Initial No. of Sperms
O0 = 0	 # Initial No. of Oocytes
E0 = 0   # Initial No. of Eggs
state = [S0,O0,E0]

if(args.model == '1'):
	## Parameters Model1  
	kg = 12
	ks = 0
	ko = 0.01
	par = [kg,ks,ko]
	model = egglayingModel1
	
elif(args.model == '2'):
## Parameters Model2
	ko = 12
	kc = 0
	kf = 0.1
	ks = 0.1
	par = [ko,kc,kf,ks]
	model = egglayingModel2
	print model
	
elif(args.model == '3'):
	
	## Parameters Model3
	ko = 12
	kf = 0.01
	par = [ko,kf]
	model = egglayingModel3
	
elif(args.model == '4'):
	## Parameters Model4
	ko = 12
	kf = 0.1
	ks = 0.1
	par = [ko,kf,ks]
	model = egglayingModel4

## Time Gird
t = np.arange(0,96,1)
## Solve the ODEs
out = odeint(model, state, t,args = (par,))
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


