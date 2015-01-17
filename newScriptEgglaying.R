##################################### Load the Libraries ###########################################
library(deSolve)
library(minpack.lm)
library(qtl)

##################################### Define necessary globals #########################################
e_times = c(2,5,27,51,72)
num_strains = 96
total_time = 96 # Hours
delta_t = 1 # Hours (for ODE)

##################################### Egg-laying models        #########################################
## Egg-laying models are called by the ode functions. Given a list of times(t), parameters(parVar),   ##
## and initial states(stateVar), these functions will return a list of State variables at each of the ##
## times.                                                                                             ##
########################################################################################################

#### Model 1  ####

egglayingModel1 = function(t,stateVar,parVar) 
        ## dE = ko*O*S
        ## dS = -dE
        ## dO = kg,i - ks*O - dE
        ## parVar = [[kg,i],ks,ko]
{
        ks = parVar[num_strains+1]
        ko = parVar[num_strains+2]
        S = c()
        with(as.list(stateVar),
{
        for(i in 1:num_strains)
        {
                
                # Find strain specific kg value
                kg = parVar[i]
                
                # Calculate dE, dS, and dO
                assign(paste("dE",i,sep=""),ko*get(paste("O",i,sep=""))*get(paste("S",i,sep="")))
                assign(paste("dS",i,sep=""),-1*get(paste("dE",i,sep="")))
                assign(paste("dO",i,sep=""),kg - ks*get(paste("O",i,sep="")) - get(paste("dE",i,sep="")))
                # Send values to output vector
                S = append(S,get(paste("dS",i,sep=""))) 
                S = append(S,get(paste("dO",i,sep="")))
                S = append(S,get(paste("dE",i,sep="")))
        }
        return(list(S))
})
}

#### Model 2  ####

egglayingModel2 = function(t,stateVar,parVar)
        ## dE = min(kf*O,ks*S)
        ## dS = -dE
        ## dO = ko,i - kc*O - dE
        ## parVar = [[kg,i],kc,kf,ks]
{
        kc = parVar[num_strains+1]
        kf = parVar[num_strains+2]
        ks = parVar[num_strains+3]
        S = c()
        with(as.list(stateVar),
{
        for(i in 1:num_strains)
        {
                # Find strain specific ko value
                ko = parVar[i]
                # Calculate dE, dS, and dO
                assign(paste("dE",i,sep=""),min(kf*get(paste("O",i,sep="")),ks*get(paste("S",i,sep=""))))
                assign(paste("dS",i,sep=""),-1*get(paste("dE",i,sep="")))
                assign(paste("dO",i,sep=""),ko - kc*get(paste("O",i,sep="")) - get(paste("dE",i,sep="")))
                # Send values to output vector
                S = append(S,get(paste("dS",i,sep=""))) 
                S = append(S,get(paste("dO",i,sep="")))
                S = append(S,get(paste("dE",i,sep="")))
        }
        return(list(S))
})
}

#### Model 3  ####

egglayingModel3 = function(t,stateVar,parVar) 
        ## dE = kf*O*S
        ## dS = -dE
        ## dO = ko,i - dE
        ## parVar = [[kg,i],ks,ko]
{
        kf = parVar[num_strains+1]
        S = c()
        with(as.list(stateVar),
{
        for(i in 1:num_strains)
        {
                
                # Find strain specific kg value
                ko = parVar[i]
                
                # Calculate dE, dS, and dO
                assign(paste("dE",i,sep=""),kf*get(paste("O",i,sep=""))*get(paste("S",i,sep="")))
                assign(paste("dS",i,sep=""),-1*get(paste("dE",i,sep="")))
                assign(paste("dO",i,sep=""),ko - get(paste("dE",i,sep="")))
                # Send values to output vector
                S = append(S,get(paste("dS",i,sep=""))) 
                S = append(S,get(paste("dO",i,sep="")))
                S = append(S,get(paste("dE",i,sep="")))
        }
        return(list(S))
})
}

#### Model 4  ####

egglayingModel4 = function(t,stateVar,parVar)
  ## dE = min(kf*O,ks*S)
  ## dS = -dE
  ## dO = ko,i - dE
  ## parVar = [[kg,i],kf,ks]
{
  kf = parVar[num_strains+1]
  ks = parVar[num_strains+2]
  S = c()
  with(as.list(stateVar),
{
  for(i in 1:num_strains)
  {
    # Find strain specific ko value
    ko = parVar[i]
    # Calculate dE, dS, and dO
    assign(paste("dE",i,sep=""),min(kf*get(paste("O",i,sep="")),ks*get(paste("S",i,sep=""))))
    assign(paste("dS",i,sep=""),-1*get(paste("dE",i,sep="")))
    assign(paste("dO",i,sep=""),ko - get(paste("dE",i,sep="")))
    # Send values to output vector
    S = append(S,get(paste("dS",i,sep=""))) 
    S = append(S,get(paste("dO",i,sep="")))
    S = append(S,get(paste("dE",i,sep="")))
  }
  return(list(S))
})
}

##################################### THE FIT FUNCTION #########################################

residualsModel = function(params, model)
        # This function is used by nls.lm to calculate the residuals for a given set of parameters. params is a list/vector
        # of parameters. Model is the function that contains the info for the egg-laying model. The parameters passed to this
        # function must be appropriate for the given model. 
{  
        #### Solve the ODE ####
        # state is a global variable containing the initial values of the state variables. Times is a global variable containing 
        # all of the times to solve for.
        out = ode(y = state, times = times, func = model, parms = params, method = solver)
        #	print(params)
        #	print(out)
        res = c()
        #### modeled - observed ####
        for(i in 1:num_strains)
        {
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (out[round(e_times/delta_t), paste("E",i,sep="")]-out[round(e_times/delta_t)-1, paste("E",i,sep="")])/delta_t
                res = append(res,modEggs-expEggs)
        }
        return (res)
}

printResult = function(params, model,modelFit)
        # This function prints the output of the least squares run
{
        #print(modelFit)
        out = ode(y = state, times = times, func = model, parms = params, method = solver)
        res = matrix(nrow=5,ncol = num_strains)
        sqres = c()
        elap_time = proc.time()-start
        print(elap_time)
        for(i in 1:num_strains)
        {	
                if(i %% 24 == 1) {
                        if(length(dev.list())!=0) {
                                dev.off()
                        }
                        file = sprintf("Plots%sto%sModel%i%s.pdf",rownames(mydata1)[i],rownames(mydata1)[min(i+23,num_strains)],choice,solver)
                        pdf(file,width = 12)
                        par(mfrow=c(4,6))                        
                }
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (out[round(e_times/delta_t), paste("E",i,sep="")]-out[round(e_times/delta_t)-1, paste("E",i,sep="")])/delta_t
                plot(times[1:length(times)-1],diff(out[,paste("E",i,sep="")]/delta_t),type="l",xlab="time",ylab="Egg-laying rate (eggs/hour)",ylim=c(0,10))
                points(e_times,expEggs)
                title(rownames(mydata1)[i])
                res[,i] = modEggs-expEggs
                sqres = append(sqres,sum(res[,i]*res[,i]))
        }
        if(length(dev.list())!=0) {
                dev.off()
        }
        outfile = sprintf("outfileModel%i%s.csv",choice,solver)
        summaryfile = sprintf("summaryModel%i%s",choice,solver)
        cat(c("Strain",rownames(mydata1),"\n"),file=outfile,sep=",")
        cat(c("ko",params[1:num_strains],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("sqres",sqres,"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res1",res[1,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res2",res[2,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res3",res[3,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res4",res[4,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res5",res[5,],"\n"),file = outfile,sep=",",append=TRUE)
        capture.output(print(summary(modelFit)), file=summaryfile)
        capture.output(print(modelFit),file=summaryfile,append=TRUE)
        capture.output(print(elap_time),file=summaryfile,append=TRUE)
        return(outfile)
}

##################################### PERFORM QTL MAPPING ##########################################
#       This function receives an argument phenoFile,which is the name of the csv file containing  #    
#       the newly generated values(ko,is and res) from the main function                           #
#       Note : You need to choose the genotype file                                                #
####################################################################################################
qtlResult = function(phenotypeFile)
{
        # Input the cross object
        print("Please Choose The Genotype File : \n")
        modelData = read.cross(format = "csvsr", genfile=file.choose(), phefile = phenoFile, genotypes = c("N2", "Het", "LSJ2"), alleles = c("N","L"))
        modelData = convert2riself(modelData)
        
        # Single Qtl scan
        modelData.scan1a = scanone(modelData, pheno.col = 2:8, method = "mr")
        modelData.perm1a = scanone(modelData, pheno.col = 2:8, n.perm = 1000, method = "mr")
        
        qtlPlotFile = sprintf("qtlPlotModel%i%s.pdf",choice,solver)
        pdf(qtlPlotFile,width=22)
        par(mfrow=c(3,7))
        
        for(i in 1:7)
        {
                plot(modelData.scan1a, lodcolumn = i, ylim = c(0,30))
                add.threshold(modelData.scan1a, perms = modelData.perm1a, alpha = 0.05, lodcolumn = i,col="red")
        }
        
        # Single Qtl scan with additive co-variate
        nurf1 = pull.geno(modelData)[,"LSJ2_N2_nurf-1"]
        modelData.scan1b = scanone(modelData,pheno.col = 2:8, method = "mr", addcovar = nurf1)
        modelData.perm1b = scanone(modelData,pheno.col = 2:8, n.perm = 1000, method = "mr", addcovar = nurf1)
        
        for(i in 1:7)
        {
                plot(modelData.scan1b, lodcolumn = i, ylim = c(0,30))
                add.threshold(modelData.scan1b, perms = modelData.perm1b, alpha = 0.05, lodcolumn = i,col="red")       
        }
        
        # Single Qtl scan with interactive co-variate
        modelData.scan1c = scanone(modelData,pheno.col = 2:8, method = "mr", addcovar = nurf1, intcov = nurf1)
        modelData.perm1c = scanone(modelData,pheno.col = 2:8, n.perm = 1000, method = "mr", addcovar = nurf1, intcovar = nurf1)
        
        for(i in 1:7)
        {
                plot(modelData.scan1c, lodcolumn = i, ylim = c(0,30))
                add.threshold(modelData.scan1c, perms = modelData.perm1c, alpha = 0.05, lodcolumn = i,col="red")    
        }        
        
        dev.off()
}

##################################### RUN MODEL #########################################
# This is a wrapper to run the fit function and print out relevant information

###################################### TAKE INPUT FROM THE USER ################################
choice = readline(prompt = "Model1---> Without Feedback Term\nModel2---> With Feedback term
Model3---> Without Feedback and Carrying Capacity\nModel4--> With Feedback and Without Carrying Capcity
           \nEnter Model Choice : (1/2/3/4)\n")
choice = as.integer(choice)

solver = readline(promp = "Enter the desired solver to be used :\n(ode45/rk4/lsoda)\n")
solver = as.character(solver)
###################################### MAIN FUNCTION ############################################
# Start the clock
start = proc.time()

# Read the Experimental data file ####
mydata1 = read.csv(file.choose(),header = TRUE,row.names=1)

# Set up the State variables
# Need an individual set of (S,O,E) for each strain. They are named as (S1,O1,E1,S2,O2,E2,....S96,O96,E96).
S = rep(300,times = num_strains)
namesS = c()
for(i in 1:num_strains)
{
  namesS = append(namesS,paste("S",i,sep=""))
}  
O = rep(0,times = num_strains)
namesO = c()
for(i in 1:num_strains)
{
  namesO = append(namesO,paste("O",i,sep=""))
}  
E = rep(0,times = num_strains)
namesE = c()
for(i in 1:num_strains)
{
  namesE = append(namesE,paste("E",i,sep=""))
}  
names(S) = namesS
names(O) = namesO
names(E) = namesE

# Generates a vector containg the initial states
states = c(S,O,E)
# Getting them in the correct order
state = c()
for(i in 1:num_strains)
{
        state = append(state,states[paste("S",i,sep="")])
        state = append(state,states[paste("O",i,sep="")])
        state = append(state,states[paste("E",i,sep="")])
}

# Set up the time interval for the experiment
delta_t = 1
times = seq(0,total_time,delta_t)

# Enable printing of the iterates
#nls.lm.control(nprint = 10)

# Set up the Parameters with initial values 
if(choice == 1)
{
        # Model1 -> (kg,i = 12, ks = 0, ko = .01)
        parametersModel1 = rep(12, times = num_strains)
        parametersModel1 = append(parametersModel1,0)
        parametersModel1 = append(parametersModel1,.01)
        lower_p = rep(0,times = num_strains+3)
        upper_p = rep(100,times = num_strains+3)
        params.fittedModel1 = nls.lm(par = parametersModel1, fn = residualsModel, lower = lower_p, upper = upper_p, jac=NULL, control=c(maxiter=1000), egglayingModel1)
        parest = coef(params.fittedModel1)
        phenoFile = printResult(parest,egglayingModel1,params.fittedModel1)
        qtlResult(phenoFile)
}
if(choice == 2)
{
        # Model2 -> (ko,i = 12, kc = 0, kf = .1, ks = .1)
        parametersModel2 = rep(12, times = num_strains)
        parametersModel2 = append(parametersModel2,0)
        parametersModel2 = append(parametersModel2,0.1)
        parametersModel2 = append(parametersModel2,0.1)
        lower_p = rep(0,times = num_strains+3)
        upper_p = rep(100,times = num_strains+3)
        params.fittedModel2 = nls.lm(par = parametersModel2, fn = residualsModel, lower = NULL,upper = NULL, jac=NULL, control=c(maxiter=1000), egglayingModel2)
        parest = coef(params.fittedModel2)
        phenoFile = printResult(parest,egglayingModel2,params.fittedModel2)
        qtlResult(phenoFile)
}
if(choice == 3)
{
        # Model3 -> (ko,i = 12, kf = .01)
        parametersModel3 = rep(10, times = num_strains)
        parametersModel3 = append(parametersModel3,.001)
        lower_p = rep(0,times = num_strains)
        lower_p = append(lower_p,.0001)
        upper_p = rep(100,times = num_strains)
        upper_p = append(upper_p,.1)
        params.fittedModel3 = nls.lm(par = parametersModel3, lower=lower_p, upper=upper_p, fn = residualsModel, jac=NULL, control=c(maxiter=1000), egglayingModel3)        
        parest = coef(params.fittedModel3)
        phenoFile = printResult(parest,egglayingModel3,params.fittedModel3)
        qtlResult(phenoFile)
}

if(choice == 4)
{
  # Model4 -> (ko,i = 12, kf = .1, ks = .1)
  parametersModel4 = rep(12, times = num_strains)
  parametersModel4 = append(parametersModel4,.1)
  parametersModel4 = append(parametersModel4,.1)
  lower_p = rep(0,times = num_strains)
  lower_p = append(lower_p,.0001)
  lower_p = append(lower_p,.0001)
  upper_p = rep(100,times = num_strains)
  upper_p = append(upper_p,1)
  upper_p = append(upper_p,1)
  
  params.fittedModel4 = nls.lm(par = parametersModel4, lower=lower_p, upper=upper_p, fn = residualsModel, jac=NULL, control=c(maxiter=1000), egglayingModel4)        
  parest = coef(params.fittedModel4)
  phenoFile = printResult(parest,egglayingModel4,params.fittedModel4)
  qtlResult(phenoFile)
}
