## Loading the required libraries
from lmfit import minimize, Parameters, Parameter,fit_report
from scipy.integrate import odeint
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np  
import csv
import argparse

########################################################################
##			Command Line Options  			      ## 
########################################################################
parser = argparse.ArgumentParser(description=
'This is a script to simulate and solve egglaying models in C.elegans')

parser.add_argument('-i','--input', help='Input file name for Experimental Data'
					,required=True)
parser.add_argument('-m','--model',help='Model No.(integer)(1/2/3/4)'
					, required=True)
parser.add_argument('-t','--time',help='Time Interval(delta_t) for the model'
					, required=True)
args = parser.parse_args()
print "Experimental Data File: %s" %(args.input)
print "Model No. Used: %s" % args.model 
print "Delta_t = %s" % args.time

########################################################################
##          Egglaying Models, take in globals containing state and    ##
##			the parameters to the model equations	      ##		  ##							              ##			
########################################################################

## Model1
def egglayingModel1(stateVar,t,params):
	## dE = ko*O*S
	## dS = -dE
	## dO = kg - ks*O - dE
	S = stateVar[0]
	O = stateVar[1]
	E = stateVar[2]
	kg = params['kg'].value
	ks = params['ks'].value
	ko = params['ko'].value
	
	dE = ko*O*S
	dS = -dE
	dO = kg - ks*O -dE
	
	return [dS,dO,dE]

## Model2
def egglayingModel2(stateVar,t,params):
	## dE = min(kf*O,ks*S)
    	## dS = -dE
    	## dO = ko - kc*O - dE
    	S = stateVar[0]
    	O = stateVar[1]
    	E = stateVar[2]
    	ko = params['ko'].value
    	kc = params['kc'].value
    	kf = params['kf'].value
    	ks = params['ks'].value
    
    	dE = min(kf*O,ks*S)
    	dS = -dE
    	dO = ko - kc*O - dE
    
    	return [dS,dO,dE]

## Model3
def egglayingModel3(stateVar,t,params):
    	## dE = kf*O*S
    	## dS = -dE
    	## dO = ko - dE
	S = stateVar[0]
	O = stateVar[1]
	E = stateVar[2]
	ko = params['ko'].value
	kf = params['kf'].value
	
	dE = kf*O*S
	dS = -dE
	dO = ko - dE
	
	return [dS,dO,dE]	
   
## Model4
def egglayingModel4(stateVar,t,params):
	## dE = min(kf*O,ks*S)
	## dS = -dE
	## dO = ko - dE
	S = stateVar[0]
	O = stateVar[1]
	E = stateVar[2]
	ko = params['ko'].value
	kf = params['kf'].value
	ks = params['ks'].value
	
	dE = min(kf*O,ks*S)
	dS = -dE
	dO = ko - dE
	
	return [dS,dO,dE]
	
########################################################################
##			     FIT FUNCTION 			      ##
##        The ODE models are called from here,returns residuals	      ##
########################################################################
def residual(params, mod, data):
	out = odeint(mod, state, t,args = (params,))
	E = out[:,2]
	eTimes = np.asarray(expTimes)
	index = np.ndarray(shape = (5,),dtype = int)
	np.around(eTimes/delta_t,out = index)
	mod = (E[index] - E[index-1])/delta_t
	modEggRate = np.array(mod,dtype = float)
	res = expPoints[i,:] - modEggRate
	return res

########################################################################
##			     PLOT FUNCTION			      ##
########################################################################
def plotEgglaying(output,mod):
	
	if(args.model == '1'):
	## Model 1 Plot : dE/dt = ko*O*S
		S = output[:,0]
		O = output[:,1]
		eggRate = ko*O*S

	elif(args.model == '2' or args.model == '4'):
	## Model 2/4 Plot : dE/dt = min(kf*O,ks*S)
		S = output[:,0]
		O = output[:,1]
		eggRate = np.minimum(kf*O,ks*S)
	
	elif(args.model == '3'):
	## Model 3 Plot : dE/dt = kf*O*S
		S = output[:,0]
		O = output[:,1]
		eggRate = kf*O*S

	## Plot the Egglaying Rate
	plt.subplot(subName)
	plt.plot(t, eggRate)
	plt.plot(expTimes, expPoints[i,0:5], marker = "o",linestyle = "--")
	plt.xlabel('Time(in hrs)')
	plt.ylabel('dE/dt')
	ptitle = "%s" %(expStrains[i])
	plt.title(ptitle)
	plt.legend(loc=0)
	
	## Residuals
	eTimes = np.asarray(expTimes)
	index = np.ndarray(shape = (5,),dtype = int)
	np.around(eTimes/delta_t,out = index)
	modRes = eggRate[index] - expPoints[i,0:5]
	modSqres = np.sum(modRes*modRes)
	print >> resFile, expStrains[i],"\t",modRes,"\t", modSqres
	

########################################################################
##			     MAIN FUNCTION                            ##
########################################################################

## Experimental File Input
f = open('/home/raghav/git/Egglaying/DataSets/Egglaying_Model.csv','rU')
csv_f = csv.reader(f)
expTimes = (2,5,27,51,72)


## Parameter and State Initialization

## States
S0 = 300 # Initial No. of Sperms
O0 = 0	 # Initial No. of Oocytes
E0 = 0   # Initial No. of Eggs
state = [S0,O0,E0]
par= Parameters()
if(args.model == '1'):
	## Parameters Model1  
	par.add('kg',value = 12)
	par.add('ks',value = 0,min = 0)
	par.add('ko',value = 0.01)
	model = egglayingModel1
	
elif(args.model == '2'):
	## Parameters Model2
	par.add('ko',value = 12,max = 100)
	par.add('kc',value = 0,min = 0)
	par.add('kf',value = 0.1)
	par.add('ks',value = 0.1,max = 2)
	model = egglayingModel2
	
elif(args.model == '3'):
	## Parameters Model3
	par.add('ko',value = 12)
	par.add('kf',value = 0.01)
	model = egglayingModel3
	
elif(args.model == '4'):
	## Parameters Model4
	par.add('ko',value = 12)
	par.add('kf',value = 0.1)
	par.add('ks',value = 0.1)
	model = egglayingModel4

## Time Gird
delta_t = float(args.time) 
t = np.arange(0,96,delta_t)

i = 0

## File handles for the summary and parameter files
#out1 = open("Report_Model%s.txt","a") %args.model
parFilename = "Parameters_Model%s.txt" % args.model
parOut = open(parFilename,"a")
expPoints = np.zeros(shape=(97,5))
expStrains = np.zeros(shape=(97),dtype="S10")
for row in csv_f:
	## Skipping first row as it has the headers
	if(i > 0):
		## Skipping first column as it has strain information
		e = np.array(row[1:6])
		expPoints[i,0:5] = e.astype(np.float)
		## Storing first column(strain information) 
		expStrains[i] = row[0]
		print expStrains[i]	
		## Call the lmfit's minimize function
		val = minimize(residual, par, args=(model,expPoints))
		#print >> out1, report_fit(par)
		## Write the fitted parameters to a file 
		if(args.model == '1'):
			print >> parOut, par['kg'].value,"\t",par['ks'].value,"\t",par['ko'].value
		elif(args.model == '2'):
			print >> parOut, par['ko'].value,"\t",par['kc'].value,"\t",par['kf'].value,"\t",par['ks'].value
		elif(args.model == '3'):
			print >> parOut, par['ko'].value,"\t",par['kf'].value
		elif(args.model == '4'):
			print >> parOut, par['ko'].value,"\t",par['kf'].value,"\t",par['ks'].value
	i = i + 1

#out1.close()
parOut.close()

fitPar = np.loadtxt(parFilename,delimiter=" ")
parFitFilename = "Parameters_Model%s_Fitted.txt" % args.model
reportFilename = "Parameters_Model%s_Report.txt" % args.model
parOut = open(parFitFilename,"a")
report = open(reportFilename,"a")

## Fitting the parameters again, some are taken as constants(after averaging them) 

## Average out and keep some parameters constant
if(args.model == '1'):
	ks = np.average(fitPar[:,1])
	ko = np.average(fitPar[:,2])

if(args.model == '2'):
	kf = np.average(fitPar[:,2])
	ks = np.average(fitPar[:,3])

if(args.model == '3'):
	kf = np.average(fitPar[:,1])

if(args.model == '4'):
	kf = np.average(fitPar[:,1])
	ks = np.average(fitPar[:,2])

## Do the fit again with the new parameter set
for i in range(1,97):
	parFit = Parameters()	
	## Write the fitted parameters to a file 
	if(args.model == '1'):	
		parFit.add('kg',value = fitPar[i-1,0])
		parFit.add('ks',value = ks,vary = False)
		parFit.add('ko',value = ko,vary = False) 
		## Call the lmfit's minimize function
		val = minimize(residual, parFit, args=(model,expPoints))
		#print >> out1, report_fit(par)
		print >> parOut, parFit['kg'].value,"\t",parFit['ks'].value,"\t",parFit['ko'].value
		print >> report, fit_report(par)
		
	elif(args.model == '2'):
		parFit.add('ko',value = fitPar[i-1,0])
		parFit.add('ks',value = ks,vary = True)
		parFit.add('kf',value = kf,vary = True)
		parFit.add('kc',value = fitPar[i-1,1])
		## Call the lmfit's minimize function
		val = minimize(residual, parFit, args=(model,expPoints))
		print >> parOut, par['ko'].value,"\t",par['kc'].value,"\t",par['kf'].value,"\t",par['ks'].value

	elif(args.model == '3'):
		parFit.add('ko',value = fitPar[i-1,0])
		parFit.add('kf',value = kf,vary = False)
		## Call the lmfit's minimize function
		val = minimize(residual, parFit, args=(model,expPoints))
		print >> parOut, par['ko'].value,"\t",par['kf'].value

	elif(args.model == '4'):
		## Call the lmfit's minimize function
		parFit.add('ko',value = fitPar[i-1,0])
		parFit.add('kf',value = kf,vary = False)
		parFit.add('ks',value = ks,vary = False)
		val = minimize(residual, parFit, args=(model,expPoints))
		print >> parOut, par['ko'].value,"\t",par['kf'].value,"\t",par['ks'].value


report.close()
parOut.close()
i = 1
j = 1
k = 1
resFilename = "Parameters_Model%s_Residuals.txt" % args.model
resFile = open(resFilename,"a")
## Read the parameters file which was created earlier 
with open(parFitFilename) as parFitFile:
	for line in parFitFile:
		if(args.model == '1'):
			kg,ks,ko = line.rstrip().split(' ')
			kg = float(kg)
			ks = float(ks)
			ko = float(ko)
			#print kg,ks,ko
			parNew = Parameters()
			parNew.add('kg',value = kg)
			parNew.add('ks',value = ks)
			parNew.add('ko',value = ko)
			outEgg = odeint(model, state, t,args = (parNew,))


		elif(args.model == '2'):
			ko,kc,kf,ks = line.rstrip().split(' ')			
			ko = float(ko)
			kc = float(kc)
			kf = float(kf)
			ks = float(ks)
			parNew = Parameters()
			parNew.add('ko',value = ko)
			parNew.add('kc',value = kc)
			parNew.add('kf',value = kf)
			parNew.add('ks',value = ks)
			outEgg = odeint(model, state, t,args = (parNew,))


		elif(args.model == '3'):
			ko,kf = line.rstrip().split(' ')
			ko = float(ko)
			kf = float(kf)
			parNew = Parameters()
			parNew.add('ko',value = ko)
			parNew.add('kf',value = kf)
			outEgg = odeint(model, state, t,args = (parNew,))


		elif(args.model == '4'):
			ko,kf,ks = line.rstrip().split(' ')
			ko = float(ko)
			kf = float(kf)
			ks = float(ks)
			parNew = Parameters()
			parNew.add('ko',value = ko)
			parNew.add('kf',value = kf)
			parNew.add('ks',value = ks)
			outEgg = odeint(model, state, t,args = (parNew,))	

		### Call the Plot Function; plot and subplot numbering,file naming and saving occurs globally 
		if(i%8 == 1):
			plt.figure(j,figsize=(15,10))
			fname = "Plot_Model%s_Strain_%stoStrain_%s" %(args.model,expStrains[i],expStrains[i+7])
			j = j + 1
			k = 1
		

		subName = "24%i" %(k)
		subName = int(subName)
		plotEgglaying(outEgg,args.model)
		if(i%8 == 0):
			plt.savefig(fname)
		i = i + 1
		k = k + 1

resFile.close()

