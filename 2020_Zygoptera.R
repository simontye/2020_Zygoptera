###############################################################
# Zygoptera phylogeny compared to Odonate Phenotypic database
# 2020.03.11
# SPT
###############################################################

# Base code from Ilvonen and Suhonen. 2016. Royal Society Open Science.
# https://royalsocietypublishing.org/doi/10.1098/rsos.160421#RSOS160421C28
# Searches for phylogenetic signals in phenotypic traits of using Pagel's λ

# Set working directory (with the tree file and datamatrix)
setwd("/Users/simontye/Documents/Research/Projects/Zygoptera/2020_Zygoptera")

# Load packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(phylolm)
library(nlme)
library(caper)

# Load phylogenetic tree
tree <- read.tree("odonataphylogeny.txt")

# Load data
data <- read.csv("speciesdata.txt", header = TRUE, sep = ";", dec = ",")
data2 <- read.csv("comparativespeciesdata.txt", header = TRUE, sep = ";", dec = ",")

###############################################################
### XXX
###############################################################

# Phylogenetical signal of weight
weight <- data[,1]
names(weight) <- rownames(data)
phylosig(tree, weight, method = "lambda", test = TRUE, nsim = 999)

# Phylogenetical signal of encapsulation
encapsulation <- data[,2]
names(encapsulation) <- rownames(data)
phylosig(tree, encapsulation, method = "lambda", test = TRUE, nsim = 999)

# Phylogenetical signal of gregarines
gregarines <- data[,3]
names(gregarines) <- rownames(data)
phylosig(tree, gregarines, method = "lambda", test=TRUE, nsim=999)

# Phylogenetical signal of water mites
watermites <- data[,4]
names(watermites) <- rownames(data)
phylosig(tree, watermites, method = "lambda", test = TRUE, nsim = 999)



###############################################################
### Phylogenetic generalized linear models of 22 different odonate species and their variables
###############################################################

#Comparative dataset
cdat <- comparative.data(data = data2, phy = tree, names.col = "Species", vcv = TRUE)

# Creating different datasets
weight        <- data2$Weight
encapsulation <- data2$Encapsulation
gregarines    <- data2$Gregarines
watermites    <- data2$Watermites

# PGLS models between two different variables
logweightencapsulation.pgls  <- pgls(log(weight) ~ encapsulation, cdat)
logweightgregarines.pgls     <- pgls(log(weight) ~ gregarines, cdat)
logweightwatermites.pgls     <- pgls(log(weight) ~ watermites, cdat)
encapsulationgregarines.pgls <- pgls(encapsulation ~ gregarines, cdat) 
encapsulationwatermites.pgls <-pgls(encapsulation~watermites, cdat)
gregarineswatermites.pgls    <- pgls(gregarines~watermites, cdat)

# Summary of the different PGLS models
summary(logweightencapsulation.pgls)
summary(logweightgregarines.pgls)
summary(logweightwatermites.pgls)
summary(encapsulationgregarines.pgls)
summary(encapsulationwatermites.pgls)
summary(gregarineswatermites.pgls)
