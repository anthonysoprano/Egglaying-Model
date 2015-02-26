library(qtl)
### Script for the qtl analysis ###
### Phenotype file used is the output csv file generated from the fit_Egglaying.py script ### 
# Input the cross object

#Get data files
print("Choose your PHENOTYPE FILE FIRST : ")
phenotypeFile = file.choose()
print("Now choose the GENOTYPE FILE :")
genotypeFile = file.choose()

#Create cross object
modelData = read.cross(format = "csvs", genfile = genotypeFile, phefile = phenotypeFile, genotypes = c("N2", "Het", "LSJ2"), alleles = c("N","L"))

#Convert into RIL
modelData = convert2riself(modelData)

#Find number of phenotypes
n_phenotypes = length(phenames(modelData))

print(n_phenotypes)

#Pull out the nurf-1 genotype for covariate analysis
nurf1 = pull.geno(modelData)[,"LSJ2_N2_nurf-1"]

#Single Qtl scan and bootstrap
modelData.scan1a = scanone(modelData, pheno.col = 2:n_phenotypes, method = "mr")
modelData.perm1a = scanone(modelData, pheno.col = 2:n_phenotypes, n.perm = 1000, method = "mr")

#QTL scan with nurf-1 as covariate
modelData.scan1b = scanone(modelData,pheno.col = 2:n_phenotypes, method = "mr", addcovar = nurf1)
modelData.perm1b = scanone(modelData,pheno.col = 2:n_phenotypes, n.perm = 1000, method = "mr", addcovar = nurf1)

# Single Qtl scan with interactive co-variate
modelData.scan1c = scanone(modelData,pheno.col = 2:n_phenotypes, method = "mr", addcovar = nurf1, intcov = nurf1)
modelData.perm1c = scanone(modelData,pheno.col = 2:n_phenotypes, n.perm = 1000, method = "mr", addcovar = nurf1, intcovar = nurf1)

qtlPlotFile = paste(dirname(phenotypeFile),"/QTLModel.pdf",sep = "")
pdf(qtlPlotFile,width=22)
par(mfrow=c(4,4))

for(i in 1:(n_phenotypes-1))
{
	print(i)
        plot(modelData.scan1a, lodcolumn = i)
        add.threshold(modelData.scan1a, perms = modelData.perm1a, alpha = 0.05, lodcolumn = i,col="black")
}
for(i in 1:(n_phenotypes-1))
{
	print(i)
	plot(modelData.scan1b, modelData.scan1c, lodcolumn = i, gap = 5, ylim = c(0,4))
        add.threshold(modelData.scan1b, perms = modelData.perm1b, alpha = 0.05, lodcolumn = i,col="black")
        add.threshold(modelData.scan1c, perms = modelData.perm1c, alpha = 0.05, lodcolumn = i,col="blue") 
}

dev.off()


