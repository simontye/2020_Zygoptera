###############################################################
# Zygoptera phylogeny compared to Odonate Phenotypic database
# 2020.03.11
# SPT
###############################################################

# Base code from Ilvonen and Suhonen. 2016. Royal Society Open Science.
# https://royalsocietypublishing.org/doi/10.1098/rsos.160421#RSOS160421C28

###############################################################

# Clear memory
rm(list=ls())

# Set working directory
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
library(phylotaR)
library(ctv)

# Update major phylogenetics packages
update.views('Phylogenetics')
#devtools::install_github('dombennett/treeman')
#devtools::install_github(repo='ropensci/phylotaR', build_vignettes=TRUE)

# Load data
tree    <- read.tree("coen_coi.tree")
opdb    <- read.csv("opdb.csv", header = TRUE, sep = ",", dec = ",")
compare <- read.csv("comparativespeciesdata.txt", header = TRUE, sep = ";", dec = ",")

###############################################################
### Create odonate phylogeny
###############################################################

wd      <- 'PATH TO WD'
ncbi_dr <- 'PATH TO NCBI BLAST TOOLS'
txid    <- 9443

# Set up working directory
setUp(wd = wd, txid = txid, ncbi_dr = ncbi_dr)

# Run pipeline
run(wd = wd)

###############################################################
### Search for phylogenetic signals in phenotypic traits using Pagel's λ
###############################################################

# Phylogenetical signal of A
A <- opdb[,1]
names(A) <- rownames(opdb)
phylosig(tree, A, method = "lambda", test = TRUE, nsim = 999)

# Phylogenetical signal of B
B <- opdb[,2]
names(B) <- rownames(opdb)
phylosig(tree, B, method = "lambda", test = TRUE, nsim = 999)

# Phylogenetical signal of C
C <- opdb[,3]
names(C) <- rownames(opdb)
phylosig(tree, C, method = "lambda", test=TRUE, nsim=999)

# Phylogenetical signal of D
D <- opdb[,4]
names(D) <- rownames(opdb)
phylosig(tree, D, method = "lambda", test = TRUE, nsim = 999)

###############################################################
### Phylogenetic generalized linear models of X odonate species and their traits
###############################################################

#Comparative dataset
cdat <- comparative.opdb(opdb = compare, phy = tree, names.col = "Species", vcv = TRUE)

# Creating different datasets
weight        <- compare$Weight
encapsulation <- compare$Encapsulation
gregarines    <- compare$Gregarines
watermites    <- compare$Watermites

# PGLS models between two different variables
logweightencapsulation.pgls  <- pgls(log(weight) ~ encapsulation, cdat)
logweightgregarines.pgls     <- pgls(log(weight) ~ gregarines, cdat)
logweightwatermites.pgls     <- pgls(log(weight) ~ watermites, cdat)

encapsulationgregarines.pgls <- pgls(encapsulation ~ gregarines, cdat) 
encapsulationwatermites.pgls <- pgls(encapsulation ~ watermites, cdat)
gregarineswatermites.pgls    <- pgls(gregarines ~ watermites, cdat)

# Summary of the different PGLS models
summary(logweightencapsulation.pgls)
summary(logweightgregarines.pgls)
summary(logweightwatermites.pgls)

summary(encapsulationgregarines.pgls)
summary(encapsulationwatermites.pgls)
summary(gregarineswatermites.pgls)
