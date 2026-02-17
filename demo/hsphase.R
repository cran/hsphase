# example of using hsphase functions
# example files simulated using HAPTRACE package
# to install the dependencies
# install.packages("gtools") 
# install.packages("snowfall")
# install.packages("Rcpp") 
# install.packages("RcppArmadillo")
# install.packages("hsphase") in case the package is not installed
# library(HAPTRACE) this library simulate genotypes and pedigree and will be uploaded to the CRAN soon
# set.seed(1)
# MAPfile = generateMAP(nChr = 2, n_markers = 1000, len_Chr = 100)
# 
# QTLfile = QTLeffects(n_QTL_Chr = 2, nChr = 2, len_Chr = 100, effect_distribution = "gamma")
# MAP_QTL_Pop = mapQTLfile(map = MAPfile, QTL = QTLfile)
# Effect = as.vector(as.numeric(MAP_QTL_Pop$Effect))
# 
# Sim_HP = generateHP(n_generations = 5, n_animals = 100, n_markers =2004, map = MAP_QTL_Pop, nSire = 20, 
#                     nDam = 50, nProgeny = 100, mutationRate = 2.5 * 10^-5, SelType = "random", Effect = Effect, 
#                     nChr = 10, h2 = 0.3, trait_mean = 40, VarE = 0.6,recL = 5, minLength = 0)
# Population_A = Sim_HP[1:(nrow(Sim_HP)/2),]
# Population_B = Sim_HP[(nrow(Sim_HP)/2 + 1):nrow(Sim_HP), ]
# 
# EffectA = Rescale(Population_A, Effect = Effect, Phenosd = 3, h2 = 0.3)
# 
# # Assign Population Code which will be reflected later in Animal IDs
# Population_A = PopulationID(Population_A, "A")
# Population_B = PopulationID(Population_B, "B")
# 
# # Recode Haplotype of Population A and B with 1 and 2 respectively
# Population_A = Recode_Haplotype(Population_A, 1)
# Population_B = Recode_Haplotype(Population_B, 2)
# 
# ### Population A
# popA <- GenOnePopulation(ngenerations=1, map=MAP_QTL_Pop,populationSim = Population_A,
#                          nSire=2,nDam=20,mutationRate=2.5 * 10^-5, SelType = "random", h2 = 0.3, trait_mean = 10, VarE = 0.6 ,Effect = Effect, recL = 5 , nChr = 2, minLength = 0, IndStartVal=1,prefixID="A",nProgeny=100)
# #> Working on the Generations: 1
# 
# ### Population B
# popB <- GenOnePopulation(ngenerations=1, map=MAP_QTL_Pop,populationSim = Population_B,
#                          nSire=5,nDam=20,mutationRate=2.5 * 10^-5, SelType = "random", h2 = 0.3, trait_mean = 30, VarE = 0.6 ,Effect = Effect, recL = 5 , nChr = 10, minLength = 0, IndStartVal=1,prefixID="B",nProgeny=100)
# #> Working on the Generations: 1
# 
# genotype <- rbind(MakeGeno(Coded_to_Haplo(popA$Haplotype)),MakeGeno(Coded_to_Haplo(popB$Haplotype)))
# genotypes = genotype
# 
# pedigree  =  rbind(popA$parents,popB$parents)
# colnames(genotype) = MAP_QTL_Pop$ID
# table(pedigree$SireID)
# MAP_QTL_Pop = MAP_QTL_Pop[,1:3]
# colnames(MAP_QTL_Pop) = c("Name","Chr","Position")
# pedigree = pedigree[,1:2]


#######################################################################
################ Running hsphase ######################################
#######################################################################
set.seed(1)
library(hsphase) # load library

# reads in a file of genotypes and a pedigree file, then splits the data into a list (of matrices - ID x SNP) of half-sibs, one for each family
# check the example files to see the format
# genotypes: ID x SNP, with rownames = ID and colnames = SNP
# genotypes: 0 - AA, 1 - AB, 2 - BB, 9 missing
# pedigree: 1st column ID, 2nd column Sire ID. Should not have a header

data(genotypes)

data(pedigree)

# split into list of half-sib groups 
halfsib <- hss(pedigree, genotypes)

# splits the family data into chromosomes
# also a list of matrices (ID x SNP): one for each chromosome X number of families
# per family, per chromosome matrices are ordered according to base pair position in the map file
# names in the list are the sireID_chromosome - can use this to find/parse subsets of data
data(map)
halfsib <- cs(halfsib,map)

# function para runs a parallel version of the main functions - it's probably what you want in practice
blocks <- para(halfsib, cpus=2, option="bmh", type = "SOCK") # bmh - builds blocks of relationship between paternal and maternal strand from sire (i.e. which chunks each offspring inherited from the sire)
sires <- para(halfsib, cpus=2, option="ssp", type = "SOCK") # ssp - phases and imputes the sire 
phased <- para(halfsib, cpus=2, option="aio", type = "SOCK") # aio - phases half-sib families (p - paternal strand from sire and m - maternal strand from dam)

# any of these function can also be run directly on a single matrix with e.g.
singleChromSirePhased <- aio(halfsib[[1]])


#######################################################################
################ ploting with hsphase #################################
#######################################################################

# visualize block structure for sire 1, chromosome 1
imageplot(blocks[[1]])

# plot number of recombinations between SNP
# needs distances between SNP - can get from map
distance <- map$Position[which(map$Chr==1)]
# plot
rplot(halfsib[[1]],distance) # plot recombinations for family one, chromosome 1


#######################################################################
################ Diagnostic tools #####################################
#######################################################################

# generates a matrix of 0/1 for recombinations
# dimension is 1 less than the number of SNP
# if a particular region cannot be resolved - all SNPs in the region are set to 1
# can help identify problems if e.g. too many recombinations are identified
recombinations <- pm(blocks[[1]])

# matches the phased haplotypes of the offspring with the phased haplotypes of its sire
# inputs are a matrix for each of the haplotypes in 0 (A), 1 (B) and 9 (missing) format
diagnostic <- hbp(phased[[1]],sires[[1]])
imageplot(diagnostic)  
# pedigree errors can be detected if there's an individual with too many recombinations
# problems with the phased data are evident if there are excessive numbers of recombinations
# this function can be used as a diagnostic tool for other phasing algorithms


#######################################################################
################ Pedigree inference and fix pedigree error ############
#######################################################################

oh <- ohg(genotypes)  # create a matrix of opposing homozygotes
hh(oh) # heatmap plot of halfsibs without pedigree colour coded sidebars

inferredPedigree <- rpoh(genotypes[, 1:ncol(halfsib[[1]])], oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3, maxRec = 5, method = "recombinations")  # infer half-sib groups pedigree 

inferredPedigree <- pedigreeNaming(inferredPedigree, pedigree) # assigns inferred pedigree to original pedigree
hh(oh, inferredPedigree, pedigree)  # heatmap with colour coded bars for the inferred and original pedigrees



#######################################################################
################### HAplotype reconstruction check ####################
#######################################################################
system.time(haplotype <- .simulateHalfsib(10, 10000, type = "haplotype"))
gMat <- groupMatSingle(haplotype$phased, 100, 2, "haplotype")
imageplot(gMat) 
example(fixSW)











