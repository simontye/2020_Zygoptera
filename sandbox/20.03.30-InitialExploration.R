setwd("~/2020_Zygoptera/")

require(ape)
require(corHMM)

# three different tree files all with the same tips. Dunno what the difference is
t <- read.nexus("data/trees/tree.nex")
#phy <- read.nexus("data/trees/Anisoptera_tree.nex")
#phy <- read.nexus("data/trees/Zygoptera_tree.nex")
d <- read.csv("data/opdb.csv")

# create a dataset
dat <- data.frame(
  sp = as.character(paste(d$Genus, d$Species, sep = "_")),
  black = black,
  plain = plain,
  ephemeral = ephemeral
)

# make sure phy and dataset are working with the same set of tips
phy <- drop.tip(t, t$tip.label[!t$tip.label %in% dat$sp])
dat <- dat[as.character(dat$sp) %in% phy$tip.label,]

# there are many doubles in the dataset. find the polymorphism and create a new state.
PolyDat <- data.frame(sp = unique(dat$sp), black = 0, plain = 0, ephemeral = 0)
for(sp in unique(dat$sp)){
  tmpDat <- dat[dat$sp %in% sp,]
  tmpDat <- apply(tmpDat[2:4], 2, function(x) paste(levels(as.factor(unique(na.omit(x)))), collapse = "/"))
  PolyDat[PolyDat$sp == sp, 2:4] <- tmpDat
}

# which of the PolyDat have empty spaces?
PolyDat <- PolyDat[!apply(apply(PolyDat, 2, function(x) x == ""), 1, any),]
PolyPhy <- drop.tip(phy, phy$tip.label[!phy$tip.label %in% PolyDat$sp])

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
            ### ### ### model fitting ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

MK <- corHMM(phy = PolyPhy, data = PolyDat, rate.cat = 1)
