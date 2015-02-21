# Todo list
# 3. Add a function to do an SVD on the data which can be output as a QTL file for mapping
# 4. Troubleshoot model 2. Is it accurate or not?

## Loading the required libraries
from scipy.integrate import ode
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
import numpy as np  
import csv
import argparse
import math
from lmfit import minimize, Parameters, Parameter, fit_report
import os
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import copy

########################################################################
##          Egglaying Models, take in globals containing state and    ##
##			the parameters to the model equations	      ##		  
##							              ##			
########################################################################

## Model1
def egglayingModel1(t,stateVar,params):
        ## Takes in a params Parameters object, keyed as = [ko1,ko2...koN, kf, kc] where N is >=1
        ## stateVar is a list [S1,O1,E1,S2,O2,E2....SN,ON,EN]
    	## dE = kf*O*S
    	## dS = -dE
    	## dO = ko - kc*O - dE

        num_strains = len(stateVar) / 3
        diffs = []

        for i in range(0,num_strains) :
                dE = params['kf'].value*stateVar[1+3*i]*stateVar[3*i]
                dS = -dE
                dO = params['ko'].value - params['kc'].value*stateVar[1+3*i] - dE
                diffs.append([dS,dO,dE])
	return diffs


## Model2
def egglayingModel2(t,stateVar,params):
        ## Takes in a params Parameters object, keyed as = [ko1,ko2...koN, kf, kc, ks] where N is >=1
        ## States are stateVar is a list [S1,O1,E1,S2,O2,E2....SN,ON,EN]
	## dE = min(kf*O,ks*S)
	## dS = -dE
	## dO = ko - kc*O - dE

        num_strains = len(stateVar) / 3
        diffs = []

        for i in range(0,num_strains) :
                dE = min(params['kf'].value*stateVar[1+3*i],params['ks'].value*stateVar[3*i])
                dS = -dE
                dO = params['ko'].value - params['kc'].value*stateVar[1+3*i] - dE
                diffs.append([dS,dO,dE])
	
	return diffs


########################################################################
##			     FIT FUNCTIONS			      ##
##        The ODE numerical differentiators and analytical functions  ##
##        are called from here, returns residuals                     ## 
########################################################################
def residual_numerical(params, model_function, data, integrator, delta_t,flag=0):
        # params is a Parameters object, 1 ko per strain. 1 kf, kc, and ks for all strains
        # model_function is a function pointer for egg-laying model
        # data is a dictionary of keys and lists [5 data points] with data['times'] = to the time points of the data
        # inital_states are [S1,O1,E1,S2,O2,E2....SN,ON,EN] where N is >=1
        # integtrator is a type of integrator used by ode
        # delta_t specifies the size of the times step
	# flag specifies the return value of this function, if 1 it returns the predicted egglaying rate, if 0 then the residuals

        initial_states = [params['S0'].value,0,0]

        r = ode(model_function).set_integrator(integrator, nsteps = 160000) 
        r.set_initial_value(initial_states,0).set_f_params(params)

        data_time_index = 0

        y = {}
        k = 0
        res = []
        pr_data=[]
        guess_data=[]

	pr_dEdt=[]

        while r.successful() and r.t < final_t and r.y[0] > .000001:
                r.integrate(r.t+delta_t)
                y[k] = r.y
                #y is [S1,O1,E1...Sn,On,En]
                k+=1

        while r.t<final_t:
                y[k] = [0, y[k-1][1], y[k-1][2]]
                k = k+1
                r.t = r.t + delta_t
        

        num_strains = len(initial_states)/3
        num_data_points = len(data['times'])
        j = 0
        print "New Function Call !"
        for key in data:
		print key
                if key != 'times' :
                        for i in range(0,num_data_points):
                                time = data['times'][i]
                                index = int(time/delta_t)
                                #Assumes that E is third data point
                                e_prime = (y[index][j+2] - y[index-1][j+2])/delta_t
                                pr_data.extend([e_prime])
                                guess_data.extend([y[index][j]*y[index][j+1]*params['kf'].value])
                                res.extend([e_prime - data[key][i]])
                        j+=1

	j = 0
	if(flag == 1):
	## Return the dE/dt values when called by the plot function
		time = np.arange(0,96,delta_t)
		for i in range(0,len(time)):		
			index = int(time[i]/delta_t)
			#Assumes that E is third data point
			dEdt = y[index][j]*y[index][j+1]*params['kf'].value
			pr_dEdt.extend([dEdt])

#        print guess_data
	if(flag==1):
		return pr_dEdt
	else:
		return res


def residual_analytic(params, data):
        # params is a Parameters object, 1 ko per strain. 1 kf, kc, and ks for all strains
        # model_function is a function pointer for egg-laying model
        # data is a dictionary of keys and lists [5 data points] with data['times'] = to the time points of the data
        # inital_states are [S1,O1,E1,S2,O2,E2....SN,ON,EN] where N is >=1
        # integtrator is a type of integrator used by ode
        # delta_t specifies the size of the times step



        ko = float(params['ko'].value)
        kf = float(params['kf'].value)
        S0 = float(params['S0'].value)
        C2 = 1.0/S0 - pow(kf/2/ko*math.pi,.5) * math.exp(kf/2/ko*pow(S0,2)) * math.erf(-1*S0*pow(kf/2/ko,.5))

       
        model = []
        all_data = []
        num_data_points = len(data['times'])
        
        for key in data:
                if key != 'times' :
                        for i in range(0,num_data_points):
                                time = data['times'][i]
                                S = (math.exp(-1*kf*(ko*pow(time,2)/2-S0*time))) /(C2 + pow(kf*math.pi/2/ko,.5) * math.exp(kf*pow(S0,2)/2/ko) * math.erf(pow(kf/2/ko,.5)*(ko*time-S0)))
                                O = ko*time - S0 + (math.exp(-1*kf*(ko*pow(time,2)/2-S0*time))) /(C2 + pow(kf*math.pi/2/ko,.5) * math.exp(kf*pow(S0,2)/2/ko) * math.erf(pow(kf/2/ko,.5)*(ko*time-S0)))
                                model.extend([kf*S*O])
                                all_data.extend([data[key][i]])

	return np.subtract(model,all_data)



########################################################################
##			  PARAMETER FUNCTION			      ##
########################################################################
def setParameters(default_params = {'ko': [12, True, 2, 50], 'kf': [.0001, True, .00000001, .01], 'kc': [0, True], 'ks': [.0001,True], 'S0':[300, True]}) : 
        # [ko, kf, kc, ks]
        parms = Parameters()
        for key in default_params:
                if len(default_params[key]) == 2:
                        parms.add(key, value = default_params[key][0], vary = default_params[key][1])
                if len(default_params[key]) == 4:
                        parms.add(key, value = default_params[key][0], vary = default_params[key][1], min = default_params[key][2], max = default_params[key][3])

        return parms

########################################################################
##			  STATE FUNCTION			      ##
########################################################################
def setStates(num_strains, default_values = [300,0,0]) : 
        states = []
        for i in range(0,num_strains) :
        	states.extend(default_values)
        return states

########################################################################
##		      INDEPENDENT FIT FUNCTION			      ##
## This function calculates the fit for all 96 strains independently  ##
########################################################################
def independent_fit_analytic(data, default_parms, int_method) :

                                                                                          
        out_data = {}
        for key in data:
                if key != 'times':
                        parms = setParameters(default_parms)
                        out_data[key] = minimize(residual_analytic, parms, method = int_method, args=({key:data[key], 'times':data['times']},))
                        if np.sum(out_data[key].residual*out_data[key].residual) > 10:
                                print key
                                #print out_data[key].params
                                #surface_plot({key:data[key], 'times':data['times']}, 4.0, 50.0, .000001, .001)


        
        return out_data

########################################################################
##		      INDEPENDENT FIT FUNCTION			      ##
## This function calculates the fit for all 96 strains independently  ##
########################################################################
def independent_fit(model, data, default_parms, int_method, integrator, delta_t) :
        out_data = {}
	# Specifying that the residual_numerical function should return the residuals, if called with flag=1 it returns the predicted egg laying rate.
	flag = 0
        for key in data:
                if key != 'times':
                        print key
                        parms = setParameters(default_parms)
                        out_data[key] = minimize(residual_numerical, parms, method = int_method, args=(model, {key:data[key], 'times':data['times']}, integrator, delta_t,flag))
        return out_data


########################################################################
##		      INDEPENDENT FIT FUNCTION			      ##
## This function   ##
########################################################################
def fix_values(model, data, default_parms, fix_parms, int_method, integrator, delta_t):
        
        new_parms = {}
        for key in default_parms:
                new_parms[key] = copy.deepcopy(default_parms[key]) 

        out_data = independent_fit(model, data, new_parms, int_method, integrator, delta_t)
        for i in fix_parms:
                parms = []
                for key in out_data:
                        parms.extend([out_data[key].params[i].value])
                new_value = np.mean(parms)
                new_parms[i][0]=new_value
                print i + " fixed with value " + str(new_value)
                new_parms[i][1]=False
                out_data = independent_fit(model, data, new_parms, int_method, integrator, delta_t)
        

        return out_data

def fix_values_analytic(data, default_parms, fix_parms, int_method):
        
        new_parms = {}
        for key in default_parms:
                new_parms[key] = copy.deepcopy(default_parms[key]) 

        out_data = independent_fit_analytic(data, new_parms, int_method)
        for i in fix_parms:
                parms = []
                for key in out_data:
                        parms.extend([out_data[key].params[i].value])
                new_value = np.mean(parms)
                new_parms[i][0]=new_value
                print i + " fixed with value " + str(new_value)
                new_parms[i][1]=False
                out_data = independent_fit_analytic(data, new_parms, int_method)
        

        return out_data

########################################################################
##			     PLOT FUNCTION			      ##
########################################################################	

def scatter_plot(obj1,obj2,directory):

        #Set up data structures
        varying_params = []
        obj1_plots = {}
        obj2_plots = {}

        #Residuals will need to be calculated and plotted
        obj1_plots['res'] = []
        obj2_plots['res'] = []

        #Create directory if needed
	if not os.path.exists(directory):
        	os.makedirs(directory)


        #Find parameters that are varying. 
        first_key = obj1.keys()[0]
        for param_key in obj1[first_key].params:
                if obj1[first_key].params[param_key].vary == True:
                        if obj2[first_key].params[param_key].vary == False:
                                print param_key + " varies in obj1 but not obj2"
                        varying_params.extend([param_key])   
                        obj1_plots[param_key] = []
                        obj2_plots[param_key] = []
        
        # Collect data that will be plotted
        for key in obj1:
                for param_key in obj1_plots:
                        if param_key == 'res':
                                #Calculate total residuals for each strain
                                obj1_plots[param_key].append(np.sum(obj1[key].residual*obj1[key].residual))
                                obj2_plots[param_key].append(np.sum(obj2[key].residual*obj2[key].residual))
                        else:
                                obj1_plots[param_key].append(obj1[key].params[param_key].value)
                                obj2_plots[param_key].append(obj2[key].params[param_key].value)

        #Create and save figures
        for param_key in obj1_plots:
                plt.figure(figsize=(15,10))
                title = "Scatter Plot parameter " + param_key
                fname = directory + "/scatter_plot_" + param_key
                plt.title(title)
                plt.scatter(obj1_plots[param_key],obj2_plots[param_key],alpha = 0.5)
                plt.savefig(fname)
                plt.close('all')
                                                        

def surface_plot(data, ko_min, ko_max, kf_min, kf_max):
        for key in data:
                if key != 'times':
                        print key
                        num_bins = 50
                        ko_step = float(ko_max - ko_min)/num_bins
                        kf_step = float(kf_max - kf_min)/num_bins
                        parms = Parameters()
                        parms.add('ko', value = ko_min, vary = True, min = 2, max = 50)
                        parms.add('kf',value = kf_min, vary = True, min = .00000001, max = .01)
                        parms.add('S0',value = 300, vary = False)
                        parms.add('kc',value = 0, vary = False)
                        parms.add('ks',value = 0, vary = False)
                        X = np.arange(ko_min,ko_max+ko_step/100,ko_step)
                        Y = np.arange(kf_min,kf_max+kf_step/100,kf_step)
                        X, Y = np.meshgrid(X, Y)
                        Z = np.zeros((len(X),len(X[0])))

                        for i in range(0,num_bins+1):
                                for j in range(0, num_bins+1):
                                        parms['ko'].set(value=X[i,j])
                                        parms['kf'].set(value=Y[i,j])
                                        res = residual_analytic(parms, {key: data[key], 'times':data['times']})
                                        Z[i,j] = sum(map(lambda i: i ** 2, res)) 
                        
                        print np.amin(np.amin(Z))
                        fig = plt.figure()
                        ax = fig.gca(projection='3d')
                        surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
                        ax.zaxis.set_major_locator(LinearLocator(10))
                        ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
                        fig.colorbar(surf, shrink=0.5, aspect=5)
                        


                        plt.show()
def derivative_analyzer(parms, data, model, integrator, d_t) :
        print "Numerical derivatives, integrator = " + integrator + " d_t = " + str(d_t)
        diff = []
        est_diff = {}
        for i in range(0,8):
                diff.extend([pow(.1,9-i)])

        for key in parms:
                if parms[key].vary:
                        for strain in data:
                                if strain != 'times':
                                        est_diff[key] = []
                                        starting_value = parms[key].value
                                        base_res = residual_numerical(parms, model, {'times':data['times'],strain:data[strain] }, integrator, d_t)
                                        base = sum(map(lambda x:x ** 2, base_res))
                                        for i in range(0,len(diff)):
                                                parms[key].set(value = starting_value + diff[i])
                                                res = residual_numerical(parms, model, {'times':data['times'],strain:data[strain]},integrator,d_t)
                                                est_diff[key].extend([base-sum(map(lambda x:x ** 2, res))])
                                        parms[key].set(value = starting_value)

        for i in range(0,len(diff)):
                out = str(diff[i])
                for key in est_diff:
                        out = out + ' ' + key + ' ' + str(est_diff[key][i]/diff[i])
                print out

def derivative_analyzer_analytic(parms, data) :
        print "Analytical derivatives"
        diff = []
        est_diff = {}
        for i in range(0,8):
                diff.extend([pow(.1,9-i)])

        for key in parms:
                if parms[key].vary:
                        for strain in data:
                                if strain != 'times':
                                        est_diff[key] = []
                                        starting_value = parms[key].value
                                        base_res = residual_analytic(parms, {'times':data['times'],strain:data[strain] })
                                        base = sum(map(lambda x:x ** 2, base_res))
                                        for i in range(0,len(diff)):
                                                parms[key].set(value = starting_value + diff[i])
                                                res = residual_analytic(parms, {'times':data['times'],strain:data[strain]})
                                                est_diff[key].extend([base-sum(map(lambda x:x ** 2, res))])
                                        parms[key].set(value = starting_value)

        for i in range(0,len(diff)):
                out = str(diff[i])
                for key in est_diff:
                        out = out + ' ' + key + ' ' + str(est_diff[key][i]/diff[i])
                print out


def large_residuals(data):
        sumsq = 0
        for key in data:
                if key != 'times':
                        mean = np.mean(data[key])
                        sumsq += sum(map(lambda x: x**2, data[key]-mean))
        print sumsq

def calc_residuals(out):
        sumsq = 0
        for key in out:
                if key != 'times':
                        sumsq += np.sum(out[key].residual*out[key].residual)
        print sumsq
	return sumsq
        		
      
def initial_analysis():
        default_params1 = {'ko': [12, True, 2, 50], 'kf': [.0004, False, .00000001, .01], 'kc': [0, False], 'ks': [.0001, False], 'S0':[300, False]}
        default_params2 = {'ko': [12, True, 2, 50], 'kf': [.0001, True, .00000001, .01], 'kc': [0, True], 'ks': [.0001, False], 'S0':[300, False]}


        print ""
        print "MODEL 1: No carrying capcity"
        default_params1 = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .01], 'kc': [0, False], 'ks': [.0001, False], 'S0':[300, False]}
        print default_params1
        out1 = independent_fit(egglayingModel1, data, default_params1, 'leastsq', 'dopri', .1)
        calc_residuals(out1)
        
        out2 = fix_values(egglayingModel1, data, default_params1, ['kf'], 'leastsq', 'dopri', .1)
        calc_residuals(out2)
        
        out3 = fix_values(egglayingModel1, data, default_params1, ['ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out3)
        
        out3 = fix_values(egglayingModel1, data, default_params1, ['kf','ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out3)
        
        print ""
        print "MODEL 1: With carrying capcity"
        default_params2 = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .01], 'kc': [0, True], 'ks': [.0001, False], 'S0':[300, False]}
        print default_params2
        out4 = independent_fit(egglayingModel1, data, default_params2, 'leastsq', 'dopri', .1)
        calc_residuals(out4)
        
        print 'kc fixed'
        out5 = fix_values(egglayingModel1, data, default_params2, ['kc'], 'leastsq', 'dopri', .1)
        calc_residuals(out5)
        
        print 'ko fixed'
        out6 = fix_values(egglayingModel1, data, default_params2, ['ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out6)
        
        print 'kf fixed'
        out7 = fix_values(egglayingModel1, data, default_params2, ['kf'], 'leastsq', 'dopri', .1)
        calc_residuals(out7)
        
        print 'kf, ko fixed'
        out8 = fix_values(egglayingModel1, data, default_params2, ['kf', 'ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out8)
        
        print 'kc, kf fixed'
        out9 = fix_values(egglayingModel1, data, default_params2, ['kc', 'kf'], 'leastsq', 'dopri', .1)
        calc_residuals(out9)
        
        print 'kc, ko fixed'
        out10 = fix_values(egglayingModel1, data, default_params2, ['kc', 'ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out10)
        
        print 'kc, kf, ko fixed'
        out11 = fix_values(egglayingModel1, data, default_params2, ['kc', 'kf','ko'], 'leastsq', 'dopri', .1)
        calc_residuals(out11)

def create_QTL_file(output, input_datafile, directory,model,initial_params,integrator,delta_t,total_residuals):
 
        varying_params = []
        output_data = {}

        #Create directory if needed
	if not os.path.exists(directory):
        	os.makedirs(directory)

        # Find parameters that are varying. 
        first_key = output.keys()[0]
        for param_key in output[first_key].params:
                if output[first_key].params[param_key].vary == True:
                        varying_params.extend([param_key])   

        output_data['ID'] = varying_params

        # Collect the data into a dictionary of rows
        for key in output:
                output_data[key] = []
                for i in range(0,len(varying_params)):
                        output_data[key].extend([output[key].params[varying_params[i]].value])
                output_data[key].extend([np.sum(output[key].residual*output[key].residual)])
                output_data[key].extend([output[key].residual[0]])
                output_data[key].extend([output[key].residual[1]])
                output_data[key].extend([output[key].residual[2]])
                output_data[key].extend([output[key].residual[3]])
                output_data[key].extend([output[key].residual[4]])

        output_data['ID'].extend(['ssq_res','res2','res5','res27','res51','res72'])

        # Open the input file to keep the strains in the right order.
        f = csv.reader(open(input_datafile,'rU'))
	
        # Open the output file and write out appropriate data
        outfile = directory + '/Modelled_phenotypes.csv'
        with open(outfile,'wb') as of:
                csvfile = csv.writer(of)
                for row in f:
                        output_data[row[0]].insert(0,row[0])
                        csvfile.writerow(output_data[row[0]])
	of.close()

	# Print the log file 
	logfile = directory + '/Modelled_phenotypes_qtllogfile'
	with open(logfile,'wb') as of:
		print >> of,"Model Used :",model,"\n"
		print >> of,"Initial Parameters :",initial_params,"\n" 
		print >> of,"Integrator Used :",integrator,"\n" 
		print >> of,"Delta_t :",delta_t,"\n"
		print >> of,"Total Residuals :",total_residuals,"\n"
	of.close()

def egglaying_plots(output,model_function,data,input_datafile,integrator,delta_t,final_time,directory):
	""" Function to plot the estimated egg-laying rate obtained with the fitted parameters, There will be a total of 96 plots, 16 in each pdf. 
	    output is the fit object returned by lmfit's minimize routine
	    model_function is a function pointer to the egg-laying model to be used
	    data is a hash keyed in with the strain vs data points and times vs experimental time points
	    input_data file is the name of the file containing the strain names and their egg-laying rate at different time points
	    integrator is the integrator we are using, viz. dopri,lsoda,etc.
	    delta_t is the time increment specified for numerical intergration
	    final_t is the final limit for the numerical integration
	    directory is the output directory for the plots
	"""
	egglaying_rate = {}
        # Create a the directory if it doesnt exist
        if not os.path.exists(directory):
		os.makedirs(directory)

	# Call the residual_numerical function with flag=1 and store the returned value of egglaying rate in a hash with keys as the strains
	# The parameters passed to this function are the ones which have been estimated by lmfit
	# The function integrates the model with these parameters
	for key in output:
		params = output[key].params
		egglaying_rate[key] = residual_numerical(params, model_function, {key:data[key], 'times':data['times']}, integrator, delta_t,flag=1)
		
	# Open the input file to keep the strains in the right order.
	f = csv.reader(open(input_datafile,'rU'))
	
	# Time grid for the plots
	t = np.arange(0,final_time,delta_t)
	
	i = 1
	j = 1
	# Plots
        for row in f:
		if(row[0]!='ID'):
			if(i%16==1):
				start_strain = row[0]
				k = 1
				plt.figure(j,figsize=(20,10))
				j = j + 1

                 	plot_title = row[0]
			plt.subplot(4,4,k)
			plt.title(plot_title)
			plt.plot(t,egglaying_rate[row[0]])
			plt.plot(data['times'], data[row[0]], marker = "o",linestyle = "--")
			plt.xlabel('Time(in hrs)')
			plt.ylabel('dE/dt')
			plt.legend(loc=0)
			if(i%16==0):
				end_strain = row[0]
				figure_name = "%s_to_%s" %(start_strain,end_strain)
				plt.savefig(directory + "/" + figure_name)
			k = k + 1
			i = i + 1
			
def calc_svd(output,input_datafile,directory):
	""" Function to calculate the Singular Value Decomposition for the residual matrix of our 96 strains. It is a 96x5 matrix
	output is the fit object returned by lmfit's minimize routine
	input_datafile is the name of the csv file containing the experimental egg-laying rates for the 96 strains
	directory is the output directory we want to save the output files to
	"""

	## SVD Theorem states : 
	## A(nxp) = U(nxn)S(nxp)Vtranspose(pxp)

	# Initialize an empty 96x5 matrix for storing the SVDs
	res_matrix = [[0 for i in range(5)] for j in range(96)]

	f = csv.reader(open(input_datafile,'rU'))
	i=0
	for row in f:
		if(row[0]!='ID'):
			res_matrix[i] = output[row[0]].residual
			i = i + 1
	
        #Create directory if needed
	if not os.path.exists(directory):
        	os.makedirs(directory)
	
	# Calculate the SVDs for the residual matrix
	U, S, V = np.linalg.svd(res_matrix,full_matrices=True)

        # Open the output file and write out appropriate data
        outfile = directory + '/svd_residuals.csv'
	
	print U
	print S
	print V
########################################################################
##			     MAIN FUNCTION                            ##
########################################################################
## Globals
final_t = 96
input_data = 'Egglaying_Model.csv'
## Experimental File Input
reader = csv.reader(open(input_data,'rU'))
data = {}
data_small = {}

#Read in data
for row in reader:
        if row[0] != 'ID' :
                data[row[0]] = map(float,row[1:])

#Add times
data['times'] = [2,5,27,51,72]

#Make small data for testing
data_small['times'] = data['times']
data_small['ln79-6'] = data['ln79-6']
data_small['ln70-7'] = data['ln70-7']

#Run code
#parms = setParameters(model1_params)
print ""
print "MODEL 1: With carrying capcity"
default_params2 = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .01], 'kc': [0, True], 'ks': [.0001, False], 'S0':[300, False]}
print default_params2
#out4 = independent_fit(egglayingModel1, data, default_params2, 'leastsq', 'dopri', .1)
#calc_residuals(out4)
        
print 'kf fixed'
out5 = fix_values(egglayingModel1, data, default_params2, ['kf'], 'leastsq', 'lsoda', .1)
tot_res = calc_residuals(out5)
#model, parameter info, integrator, delta_t, total residuals
create_QTL_file(out5, input_data, "./QTLData","EgglayingModel1",default_params2,'lsoda',.1,tot_res)
egglaying_plots(out5,egglayingModel1,data,input_data,'dopri',.1,96,"./ModelPlots")
calc_svd(out5,"Egglaying_Model.csv","./QTLData")

#derivative_analyzer_analytic(parms, data_small)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dop853', .01)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dop853', .1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dop853', 1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dopri', .01)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dopri', .1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'dopri', 1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'vode', .01)
#derivative_analyzer(parms, data_small, egglayingModel1, 'vode', .1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'vode', 1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'lsoda', .01)
#derivative_analyzer(parms, data_small, egglayingModel1, 'lsoda', .1)
#derivative_analyzer(parms, data_small, egglayingModel1, 'lsoda', 1)



"""
print "analytical"
model1_params = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .001], 'kc': [0, False], 'ks': [0, False], 'S0':[300, True, 250,500]}
print 'kf fixed'
out_a_f = fix_values_analytic(data, model1_params, ['kf'], 'leastsq')
calc_residuals(out_a_f)
create_QTL_file(out_a_f, input_data, 'Model1_Analytical_kf_fixed_ko_S0_varies')
print 'kf,S0 fixed'
#out_a_fS = fix_values_analytic(data, model1_params, ['kf','S0'], 'leastsq')
#calc_residuals(out_a_fS)
#create_QTL_file(out_a_fS, input_data, 'Model1_Analytical_kfS0_fixed_ko_varies')



print "MODEL 1: With carrying capcity"
#default_params2 = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .01], 'kc': [0, False], 'ks': [.0001, False], 'S0':[300, True]}

print 'kf fixed'
#out_kf = fix_values(egglayingModel1, data, default_params2, ['kc'], 'leastsq', 'dopri', .1)
#calc_residuals(out_kf)
#create_QTL_file(out_kf, input_data, 'Model1_dopri_kfkc_fixed_ko_S0_varies')
#print 'kf,kc fixed'
#out_kfkc = fix_values(egglayingModel1, data, default_params2, ['kf','kc'], 'leastsq', 'dopri', .1)
#calc_residuals(out_kfkc)
#create_QTL_file(out_kfkc, input_data, 'Model1_dopri_kf_kc_fixed_ko_varies')

#print "MODEL 2: No carrying capcity"
#default_params3 = {'ko': [12, True, 2, 50], 'kf': [.03, True, .00000001, .1], 'kc': [0, False], 'ks': [.1, True], 'S0':[300, False]}
#print 'kf fixed'
#out5 = fix_values(egglayingModel2, data, default_params3, ['kf'], 'leastsq', 'dopri', .1)
#calc_residuals(out5)
#create_QTL_file(out5, input_data, 'Model2_dopri_kf_fixed_ko_ks_varies')
#print 'kf,ks fixed'
#out6 = fix_values(egglayingModel2, data, default_params3, ['kf','ks'], 'leastsq', 'dopri', .1)
#calc_residuals(out6)
#create_QTL_file(out6, input_data, 'Model2_dopri_kf_ks_fixed_ko_varies')


#print "dopri .1"
#out_dopri = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'dopri', .1)
#calc_residuals(out_dopri)
#scatter_plot(out_a,out_dopri,'M1_VarykokfS0_analaticvsdopri')


print "dopri .1"
out_dopri = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'dopri', .1)
calc_residuals(out_dopri)
print "vode .1"
out_vode = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'vode', .1)
calc_residuals(out_vode)
print "lsoda .1"
out_lsoda = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'lsoda', .1)
calc_residuals(out_lsoda)

model1_params = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .001], 'kc': [0, True], 'ks': [0, False], 'S0':[300, False]}
print "kc dop853 .1"
out_dop853_kc = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'dop853', .1)
calc_residuals(out_dop853_kc)
print "kc dopri .1"
out_dopri_kc = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'dopri', .1)
calc_residuals(out_dopri_kc)
print "kc vode .1"
out_vode_kc = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'vode', .1)
calc_residuals(out_vode_kc)
print "kc lsoda .1"
out_lsoda_kc = independent_fit(egglayingModel1, data, model1_params, 'leastsq', 'lsoda', .1)
calc_residuals(out_lsoda_kc)


#default_params1 = {'ko': [12, True, 2, 50], 'kf': [.03, True, .00000001, 1], 'kc': [0, False], 'ks': [.5, True], 'S0':[300, False]}
#out1 = independent_fit(egglayingModel2, data, default_params1, 'leastsq', 'dop853', .1)
#out2 = independent_fit(egglayingModel2, data, default_params1, 'leastsq', 'dopri', .05)
#default_params1 = {'ko': [12, True, 2, 50], 'kf': [.0004, True, .00000001, .01], 'kc': [0, True], 'ks': [.0001, False], 'S0':[300, False]}
#out3 = independent_fit(egglayingModel1, data, default_params1, 'leastsq', 'dopri', .1)


#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .1)
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .05)
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .02)
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .01) 
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .005)
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .002) 
#derivative_analyzer(parms, data_small, egglayingModel2, 'dop853', .001) 




#out1 = independent_fit(egglayingModel1, data_small, [12,.0001,0,.0001], [True,False,False,False], 'leastsq', 'lsoda', .01)
#out1 = independent_fit(egglayingModel1, data, [12,.0001,0,.0001], [True,True,False,False], 'leastsq', 'dopri', .1)
#print "Residuals, no model"
#large_residuals(data)


scatter_plot(out_a,out_dop853,'M1_Nokc_analaticvsdop853')
scatter_plot(out_a,out_dopri,'M1_Nokc_analaticvsdopri')
scatter_plot(out_a,out_vode,'M1_Nokc_analaticvsvode')
scatter_plot(out_a,out_lsoda,'M1_Nokc_analaticvslsoda')

scatter_plot(out_dop853_kc,out_dopri_kc,'M1_kc_dop853vsdopri')
scatter_plot(out_dop853_kc,out_vode_kc,'M1_kc_dop853vsvode')
scatter_plot(out_dop853_kc,out_lsoda_kc,'M1_kc_dop853vslsoda')
"""
#surface_plot(data_small, 1.0, 50.0, .00001, .0005);
