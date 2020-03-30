setwd("~/2020_Zygoptera/")

require(ape)
require(corHMM)

# three different tree files all with the same tips. Dunno what the difference is
phy <- read.nexus("data/trees/tree.nex")
#phy <- read.nexus("data/trees/Anisoptera_tree.nex")
#phy <- read.nexus("data/trees/Zygoptera_tree.nex")
dat <- read.csv("data/opdb.csv")

dat