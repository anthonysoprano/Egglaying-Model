########################################################################################################
##      This script solves the egg-laying dynamics in C.elegans and fits the parameters so            ##
##      that the residual sum of squares is minimized w.r.t the experimental data                     ##
##                                                                                                    ##
##      QTL mapping is then performed using these parameters and resiudals as the phenotypes          ##
##      on the C.elegans genome                                                                       ##
########################################################################################################

# Load the required libraries 
library(deSolve)        
library(minpack.lm)
library(qtl)

######################################       Egg-laying models        ###################################
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
        return(list(c(dS,dO,dE)))
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
        return(list(c(dS,dO,dE)))
})
}

# Model 3 

egglayingModel3 = function(t,stateVar,parVar) 
        ## dE = kf*O*S
        ## dS = -dE
        ## dO = ko - dE
        ## parVar = [ko,kf]
{
        with(as.list(c(stateVar,parVar)),
{
        dE = kf*O*S
        dS = -dE
        dO = ko - dE
        return(list(c(dS,dO,dE)))
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
        return(list(c(dS,dO,dE)))
})
}

######################################### The Fit Function ################################################
residualsModel = function(params, model)
{       # This function is used by nls.lm to calculate the residuals for a given set of parameters. params is a list/vector
        # of parameters. Model is the function that contains the info for the egg-laying model. The parameters passed to this
        # function must be appropriate for the given model. 
        
        ## Solve the ODEs
        # state is a global variable containing the initial values of the state variables. Times is a global variable containing 
        # all of the times to solve for.
        out = ode(y = state, times = times, func = model, parms = params, method = solver)
        expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
        modEggs = (out[round(e_times/delta_t), "E"]-out[round(e_times/delta_t)-1, "E"])/delta_t
        res = modEggs-expEggs
        #print(expEggs)
        return(res)
}

####################################### The Plot Function ##################################################
printResult = function(params, model,modelFit)
        # This function prints the output of the least squares run
{
        #print(modelFit)
        out = ode(y = state, times = times, func = model, parms = params, method = solver)
        
        ## Plot 
        if(i%%24 == 1)
        {
                #x11()
                file = sprintf("Plot_Model_%i_Strain_%stoStrain_%s.pdf",choice,rownames(mydata1[i]),rownames(mydata1)[i+23])
                pdf(file,width = 12)
                par(mfrow=c(4,6))
        }
        plot(times[1:length(times)-1],diff(out[,"E"]/delta_t),type="l",xlab="time",ylab="Egg-laying rate (eggs/hour)",ylim=c(0,10))
        expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
        #modEggs = (out[round(e_times/delta_t), "E"]-out[round(e_times/delta_t)-1, "E"])/delta_t
        points(e_times,expEggs)
        title(rownames(mydata1)[i])
        #res = modEggs - expEggs
        if(i%%24 == 0)
        {
                dev.off()
        }
        outfile = sprintf("outfileModel%i%s.csv",choice,solver)
        summaryfile = sprintf("summaryModel%i%s",choice,solver)

        capture.output(print(summary(modelFit)), file=summaryfile,append=TRUE)
        capture.output(print(modelFit),file=summaryfile,append=TRUE)
        capture.output(print(rownames(mydata1)[i]),file = summaryfile,append=TRUE)
        #capture.output(print(elap_time),file=summaryfile,append=TRUE)
        return (out)
}
phenoFileGenerator = function()
{
        
        outfile = sprintf("outfileModel%i%s.csv",choice,solver)
        summaryfile = sprintf("summaryModel%i%s",choice,solver)
        cat(c("Strain",rownames(mydata1),"\n"),file=outfile,sep=",")
        if(choice == 1){
                cat(c("kg",params_kg[1:num_strains],"\n"),file = outfile,sep=",",append=TRUE)
        }
        else{
                cat(c("ko",params_ko[1:num_strains],"\n"),file = outfile,sep=",",append=TRUE) 
        }
        cat(c("sqres",sqres,"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res1",res[1,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res2",res[2,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res3",res[3,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res4",res[4,],"\n"),file = outfile,sep=",",append=TRUE)
        cat(c("res5",res[5,],"\n"),file = outfile,sep=",",append=TRUE)     
}

##################################### PERFORM QTL MAPPING ##########################################
##         This function receives an argument phenoFile,which is the name of the csv file         ##    
##         containing the newly generated values(ko and res) from the main function               ##
##         Note : You need to choose the genotype file                                            ##
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

#############################################  RUN MODEL ###################################################
# This is a wrapper to run the fit function and print out relevant information

###################################### TAKE INPUT FROM THE USER ############################################
choice = readline(prompt = "Model1---> Without Feedback Term\nModel2---> With Feedback term
Model3---> Without Feedback and Carrying Capacity\nModel4--> With Feedback and Without Carrying Capcity
           \nEnter Model Choice : (1/2/3/4)\n")
choice = as.integer(choice)
solver = readline(prompt = "Enter the desired solver to be used :\n(ode45/rk4/lsoda)\n")
solver = as.character(solver)
###################################### MAIN FUNCTION #######################################################
## Start the clock
start = proc.time()

## Read the Experimental data file ####
#mydata1 = read.csv(file.choose(),header = TRUE,row.names=1)

## No. of Strains
num_strains = 24
## Initialize State Variables
state = c(S=300,O=0,E=0)
## Set up the time interval for the experiment
total_time = 96
delta_t = 1
e_times = c(2,5,27,51,72)
times = seq(0,total_time,delta_t)
model = "egglayinModel%i"
# Set up the Parameters with initial values 
res = matrix(nrow=5,ncol = num_strains)
sqres = c()
if(choice == 1)
{
        # Model1 -> (kg,i = 12, ks = 0, ko = .01)
        parametersModel1 = c(kg = 12, ks = 0.01, ko = 0.01)
        ## Initialize empty vectors for storage
        params_kg = rep(0,times = num_strains)
        params_ks = rep(0,times = num_strains)
        params_ko = rep(0,times = num_strains)
        pkg = c()
        
        ## Do the first Fit
        for(i in 1:num_strains)
        {
                params.fittedModel1 = nls.lm(par = parametersModel1, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel1)
                parest = coef(params.fittedModel1)
                #print(parest)
                pkg = append(pkg,parest[kg])
                params_ks[i] = parest[[2]]
                params_ko[i] = parest[[3]]

        }
        ## Average some parameters
        ks = mean(params_ks)
        ko = mean(params_ko)
        
        ## Do the fit again, note : par is the new parameter variable, it excludes the above variables
        for(i in 1:num_strains)
        {
                parFree = pkg[i]
                params.fittedModel1 = nls.lm(par = parFree, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel1)
                parest = coef(params.fittedModel1)
                params_kg[i] = parest[[1]]
                output = printResult(parest,egglayingModel1,params.fittedModel1)
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (output[round(e_times/delta_t), "E"]-output[round(e_times/delta_t)-1, "E"])/delta_t
                res[,i] = modEggs-expEggs
                sqres = append(sqres,sum(res[,i]*res[,i]))
        }
        phenoFile = phenoFileGenerator()
        #qtlResult(phenoFile)
}
## Empty vector for storing ko values (to be outputted to file)
pko = c()
if(choice == 2)
{
        # Model2 -> (ko,i = 12, kc = 0, kf = .1, ks = .1)
        parametersModel2 = c(ko = 12, kc = 0.1, kf = 0.1, ks = 0.1)
        params_ko = rep(0,times = 96)
        params_kc = rep(0,times = 96)
        params_kf = rep(0,times = 96)
        params_ks = rep(0,times = 96)
        lower_p = rep(0,times = 4)
        upper_p = rep(100,times = 4)
        for(i in 1:num_strains)
        {
                params.fittedModel2 = nls.lm(par = parametersModel2, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel2)
                parest = coef(params.fittedModel2)
                #print(parest)
                pko = append(pko,parest[[1]])
                params_kc[i] = parest[[2]]
                params_kf[i] = parest[[3]]
                params_ks[i] = parest[[4]]

        }
        kf = mean(params_kf)
        ks = mean(params_ks)
        kc = mean(params_kc)
        for(i in 1:num_strains)
        {
                parFree = pko[i]
                params.fittedModel2 = nls.lm(par = parFree, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel2)
                parest = coef(params.fittedModel2)
                params_ko[i] = parest[[1]]
                output = printResult(parest,egglayingModel2,params.fittedModel2)
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (output[round(e_times/delta_t), "E"]-output[round(e_times/delta_t)-1, "E"])/delta_t
                res[,i] = modEggs-expEggs
                sqres = append(sqres,sum(res[,i]*res[,i]))
                
        }
        phenoFile = phenoFileGenerator()
        #qtlResult(phenoFile)
}
if(choice == 3)
{
        # Model3 -> (ko = 12, kf = .01)
        parametersModel3 = c(ko = 12, kf = 0.001)
        params_ko = rep(0,times = 96)
        params_kf = rep(0,times = 96)
        lower_p = c(0,0.0001)
        upper_p = c(100,0.1)
        for(i in 1:num_strains)
        {
                params.fittedModel3 = nls.lm(par = parametersModel3, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel3)
                parest = coef(params.fittedModel3)
                pko = append(pko,parest[[1]])
                params_kf[i] = parest[[2]]
                #print(parest)                
        }
        kf = mean(params_kf)
        for(i in 1:num_strains)
        {
                parFree = pko[i]
                params.fittedModel3 = nls.lm(par = parFree, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel3)
                parest = coef(params.fittedModel3)
                params_ko[i] = parest[[1]]
                output = printResult(parest,egglayingModel3,params.fittedModel3)
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (output[round(e_times/delta_t), "E"]-output[round(e_times/delta_t)-1, "E"])/delta_t
                res[,i] = modEggs-expEggs
                sqres = append(sqres,sum(res[,i]*res[,i]))
        }
        phenoFile = phenoFileGenerator()
        #qtlResult(phenoFile))
}

if(choice == 4)
{
        # Model4 -> (ko = 12, kf = .1, ks = .1)
        parametersModel4 = c(ko = 12, kf = 0.1, ks = 0.1)
        params_ko = rep(0,times = 96)
        params_kf = rep(0,times = 96)
        params_ks = rep(0,times = 96)
        lower_p = c(0,0.0001,0.0001)
        upper_p = c(100,1,1)
        for(i in 1:num_strains)
        {
                params.fittedModel4 = nls.lm(par = parametersModel4, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel3)
                parest = coef(params.fittedModel4)
                #print(parest)
                pko = append(pko,parest[1])
                params_kf[i] = parest[[2]]
                params_ks[i] = parest[[3]]
        }
        
        kf = mean(params_kf)
        ks = mean(params_ks)
        
        for(i in 1:num_strains)
        {
                parFree = pko[i]
                params.fittedModel4 = nls.lm(par = parFree, fn = residualsModel, lower = NULL, upper = NULL, jac = NULL, control=c(maxiter=1000), egglayingModel4)
                parest = coef(params.fittedModel4)
                params_ko[i] = parest[[1]]
                output = printResult(parest[1],egglayingModel4,params.fittedModel4)
                expEggs = c(mydata1[i,]$T1,mydata1[i,]$T2,mydata1[i,]$T3,mydata1[i,]$T4,mydata1[i,]$T5)
                modEggs = (output[round(e_times/delta_t), "E"]-output[round(e_times/delta_t)-1, "E"])/delta_t
                res[,i] = modEggs-expEggs
                sqres = append(sqres,sum(res[,i]*res[,i]))
        }
        phenoFile = phenoFileGenerator()
        #qtlResult(phenoFile)
}

