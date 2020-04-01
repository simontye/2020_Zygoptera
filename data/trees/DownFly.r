# This is how we build our tree...

# Download ----------> 

if(FALSE) { # 1. # DownFly # download (almost) every gene with every species
load(file = "C:/Users/Owner/Desktop/BodySize/DownFly/Species.rda")

library(seqinr)

package = function(sequences,annots) { # package up genetic data into lists
hit = list()

length_annots = length(annots)
if(length_annots >= 500) length_annots = 500 # check

for(i in 1:length_annots) {
Seq = sequences[[i]]
Annot = annots[[i]]
hit[[i]] = list(Annot=Annot,Seq=Seq)
}
return(hit)
}

library(seqinr)
sp = Species$GenusSpecies
SL = split(sp, ceiling(seq_along(sp)/1000)) # split into chunks to download separately

DownFly = function(sp) { # Main DownFly function

for(i in 1:length(sp)) {
choosebank("genbank", timeout = 30)
query = paste0("SP=",sp[i])
Request = try(query(listname = "genes", query = query, socket = autosocket(), invisible = FALSE, verbose = TRUE, virtual = FALSE))

if(!class(Request) == "try-error") {
reqLength = length(genes$req) # get request length
print(paste0("Num Records: ",reqLength))

if(reqLength < 500) {
sequences = getSequence(genes$req) 
annots = getAnnot(genes$req, nbl = 100)
hits = package(sequences,annots)
} else hits = "Too Many Records"


} else hits = NA

print(sp[i])
print(hits)
closebank()

save(hits,file = paste0("C:/Users/Owner/Desktop/bodysize/DownFly/hits/",sp[i],".rda")) # save
}

}

# Run DownFly in parallel

library(foreach)
library(doParallel)
library(random)

cl = makePSOCKcluster(length(SL))             
clusterSetRNGStream(cl)
registerDoParallel(cl)

X = foreach(i=1:length(SL), .packages=c('seqinr'), .multicombine=TRUE) %dopar% {  
sp = SL[[i]]
DownFly(sp)
}
stopCluster(cl)

}

if(FALSE) { # 2. # Process hits into genes.rda which is raw genes with annotations
Files = list.files("C:/Users/Owner/Desktop/BodySize/DownFly/hits/")

Names = gsub(".rda","",Files) # names for list
genes = list()
for(i in 1:length(Files)) { 
File = Files[i]
load(paste0("C:/Users/Owner/Desktop/BodySize/DownFly/hits/",File))
genes[[i]] = hits
}
names(genes) = Names
save(genes,file="C:/Users/Owner/Desktop/BodySize/DownFly/genes.rda") # save raw genes
}

if(FALSE) { # 3. # Process genes.rda into dirty_products_1.rda data.frame
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/genes.rda")

CleanProduct = function(Product) { # Clean up Annots
Product = gsub("\"","",Product) # remove quotes
Product = gsub("/","",Product) # remove slash
Product = gsub("product=|note=","",Product) # words
Product = gsub("^\\s+", "", Product) # trim leading whitespace
return(Product)
}

CleanRegion = function(Region) {
Region = gsub("[A-Z]|[a-z]","",Region) # remove letters
Region = gsub(" ","",Region) # remove spaces
Region = gsub(">","",Region) # remove
Region = gsub("<","",Region) # remove 
Region = gsub("_","",Region) # remove 
return(Region)
}


Check = function(spLevel) { # check if 
if(length(spLevel[[1]]) > 1) {
out = TRUE
} else out = FALSE
return(out)
}

GetBlocks = function(Annot) {

SourceLines = Annot[which(grepl(" source ",Annot)):length(Annot)]
Where = which(grepl("\\.\\.",SourceLines))
Blocks = list()
for(i in 1:(length(Where))) {
if(i < length(Where)) Blocks[[i]] = SourceLines[Where[i]:(Where[i+1]-1)]
else Blocks[[i]] = SourceLines[Where[length(Where)]:length(SourceLines)]
}

return(Blocks)
}

ProcessBlocks = function(Blocks) {

Regions = c()
Products = c()

for(i in 1:length(Blocks)) {
Block = Blocks[[i]] 
Product = Block[grepl("product=|note=",Block)]
Region = Block[grepl("\\.\\.",Block)]

Region = CleanRegion(Region)
Product = CleanProduct(Product)

Regions[i] = Region
if(length(Product) != 0) Products[i] = Product else Products[i] = NA

}


D = data.frame(Products,Regions)
return(D)
}


sp = names(genes)
DDList = list()
for(i in 1:length(genes)) {
spLevel = genes[[i]]
Name = sp[i]
DList = list()

if(Check(spLevel)) {

for(j in 1:length(spLevel)) {
Annot = spLevel[[j]]$Annot
Seq = paste(spLevel[[j]]$Seq,collapse = "")
Blocks = try(GetBlocks(Annot))

if(!class(Blocks) == "try-error") {
Blocks = Blocks[2:length(Blocks)]
} else Blocks = NA

if(!is.na(Blocks)) {
D = ProcessBlocks(Blocks)
DList[[j]] = cbind(Name,D,Seq)
}
}

} else DList[[1]] = data.frame(Name,Products = NA,Regions = NA, Seq = NA)
DDList[[i]] = do.call("rbind",DList)

}

dirty_products = do.call("rbind",DDList) 

save(dirty_products,file="C:/Users/Owner/Desktop/BodySize/DownFly/dirty_products_1.rda") # first stage of dirty products

}

if(FALSE) { # 4. # clean up dirty_products_1.rda into dirty_prodcuts_2.rda

load(file="C:/Users/Owner/Desktop/BodySize/DownFly/dirty_products_1.rda")
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Species.rda")

Products = dirty_products # rename 

# clean up 
Regions = Products$Regions
Regions = gsub("=","",Regions)
Regions = gsub(":","",Regions)
Regions = gsub(",,","",Regions)
Regions = gsub("\\(","",Regions)
Regions = gsub("\\)","",Regions)
Regions = gsub("\\/","",Regions)

Products$Regions = Regions
Products = na.omit(Products) # remove missing

GetSeqChunks = function(chunks,Seq) { # get chunks
SeqChunks = list()
for(i in 1:length(chunks)) {
Start = as.numeric(chunks[[i]][1])
Finish = as.numeric(chunks[[i]][2])
SeqChunks[[i]] = substring(Seq,Start,Finish)
}

return(SeqChunks)
}

L = strsplit(Products$Regions, "\\.\\.|\\,")
NewSeq = c()
for(i in 1:length(L)) {
Region = L[[i]]
chunks = split(Region, ceiling(seq_along(Region)/2)) # split into chunks of two
Seq = Products$Seq[i]
NewSeq[i] = paste(unlist(GetSeqChunks(chunks,Seq)),collapse = "")
}

nbp = nchar(NewSeq)
Products$nbp = nbp
Products$FullSeq = Products$Seq
Products$Seq = NewSeq
Products$GenusSpecies = Products$Name


Products = merge(Species,Products, by = "GenusSpecies") # add full taxonomic information

dirty_products = Products
# save(dirty_products,file="C:/Users/Owner/Desktop/BodySize/DownFly/dirty_products_2.rda") # first stage of dirty products


}

if(FALSE) { # 5. # rename product names; process into products for alignment; dirty_prodcuts_2.rda --> Products.rda; DNA with gene names
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/dirty_products_2.rda") # second stage

Products = dirty_products # rename
Products$ProductsOrginal = Products$Products # save original name before replacing synonyms

# replace names
p = Products$Products # load into vector p to process

Rsyn = function(synonyms,replacement,p) { # remove synonym function

for(i in 1:length(synonyms)) {
p = gsub(synonyms[i],replacement,p)
}
return(p)
}

# ReplaceSynonyms
replacement = "cytochrome oxidase I"
synonyms = c("cytochrome c oxidase I", "cytochrome c oxidase I subunit","cytochrome c oxidase subunit I", "cytochrome oxidase I", "cytochrome oxidase subunit 1", "cytochrome oxidase subunit I","COI","cytochrome oxidase I subunit","similar to cytochrome oxidase I")
p = Rsyn(synonyms,replacement,p)

replacement = "cytochrome oxidase II"
synonyms = c("cytochrome c oxidase subunit 2", "cytochrome oxidase II", "cytochrome oxidase subunit 2", "cytochrome oxidase subunit II") 
p = Rsyn(synonyms,replacement,p)
 
replacement = "18S ribosomal RNA"
synonyms = c("18S ribosomal RNA","18S rRNA'", "contains 18S ribosomal RNA, internal transcribed", "18S small subunit ribosomal RNA", "contains 18S ribosomal RNA, internal transcribed")
p = Rsyn(synonyms,replacement,p)

replacement = "internal transcribed spacer 1"
synonyms = c("internal transcribed spacer 1","internal transcribed spacer 1, ITS1")
p = Rsyn(synonyms,replacement,p)

replacement = "internal transcribed spacer 2"
synonyms = c("internal transcribed spacer 2","internal transcribed spacer 2, ITS2", "internal transcribed spacer 2, type 2")
p = Rsyn(synonyms,replacement,p)

replacement = "histone H3"
synonyms = c("histone 3", "histone H3")
p = Rsyn(synonyms,replacement,p)

replacement = "NADH dehydogenase subunit I"
synonyms = c("NADH dehydogenase subunit I", "NADH dehydrogenase subunit 1", "NADH dehydrogenase subunit 1; putative nuclear","NADH dehydrogenase subunit I; coding region not","ND1")
p = Rsyn(synonyms,replacement,p)

replacement = "elongation factor 1 alpha"
synonyms = c("elongation factor-1 alpha","elongation factor 1-alpha","elongation factor 1 alpha")
p = Rsyn(synonyms,replacement,p)

replacement = "myosin light chain"
synonyms = c("mlc","myosin","myosin light chain")
p = Rsyn(synonyms,replacement,p)

replacement = "16S ribosomal RNA"
synonyms = c("contains 16S ribosomal RNA and tRNA-Leu")
p = Rsyn(synonyms,replacement,p)

Products$Products = p # put p back into data.frame
# save(Products,file="C:/Users/Owner/Desktop/BodySize/DownFly/Products.rda") # second stage

}

# Alignment ------------->

if(FALSE) { # 6. # Makes Fastas for Muscle

# Align Products.rda
load(file = "C:/Users/Owner/Desktop/BodySize/DownFly/Products.rda")

# 12S 16S COI COII 18S 28S EF1 H3 # genes used by Ware
Genes = c("tRNA-Val","histone H3","12S ribosomal RNA","NADH dehydogenase subunit I","tRNA-Leu","cytochrome oxidase II","elongation factor 1 alpha","internal transcribed spacer 2","5.8S ribosomal RNA","internal transcribed spacer 1","18S ribosomal RNA","28S ribosomal RNA","16S ribosomal RNA","cytochrome oxidase I")

Products = Products[Products$Products %in% Genes,] # filter products

# Filter Outliers
PL = split(Products,Products$Products) # Products List

L = list(); for(i in 1:length(PL)) L[[i]] = PL[[i]][PL[[i]]$nbp/max(PL[[i]]$nbp) > 0.2,] # remove outliers
names(L) = names(PL)

Products = do.call("rbind",L) # combine 

# Outgroup
load(file = "C:/Users/Owner/Desktop/BodySize/DownFly/outgroup.rda")

Products = rbind(Products,outgroup) # put in outgroup

GenesList = list(); for(i in 1:length(Genes)) GenesList[[i]] = Products[Products$Products == Genes[i],]

GetSequencesFromFocal = function(FocalGene) { #

Seq = FocalGene$Seq
Names = FocalGene$Name
sequences = list()
for(i in 1:length(Names)) {

S = tolower(substring(Seq[i], seq(1,nchar(Seq[i]),1), seq(1,nchar(Seq[i]),1))) # split into ind characaters
sequences[[i]] = gsub("n","",S) # remove nnnnnn from genes
}

names(sequences) = Names

return(sequences)
}

GenesSequences = list()
for(i in 1:length(GenesList)) {
FocalGene = na.omit(GenesList[[i]]) # remove missing ends up removing outgroup
FocalGene = FocalGene[!duplicated(FocalGene$Name),] # remove duplicated species
GenesSequences[[i]] = GetSequencesFromFocal(FocalGene)
}
names(GenesSequences) = Genes

for(i in 1:length(Genes)) {
sequences = GenesSequences[[i]]
Dir = "C:/Users/Owner/Desktop/BodySize/DownFly/Fastas/"
File = paste0(Dir,Genes[i],".fasta")
# seqinr::write.fasta(sequences, names(sequences), file=File)
}

}

if(FALSE) { # 7a. # Run Muscle Alignment in Parallel
Dir = "C:/Users/Owner/Desktop/BodySize/DownFly/Fastas/"
Files = list.files(Dir)

library(foreach)
library(doParallel)
library(random)
library(muscle)

cl = makePSOCKcluster(10)             
clusterSetRNGStream(cl)
registerDoParallel(cl)

X = foreach(i=1:length(Files), .packages=c('muscle'), .multicombine=TRUE) %dopar% {  
File = paste0(Dir,Files[i])
aln = muscle(seqs = File)
# save(aln,file = paste0("C:/Users/Owner/Desktop/BodySize/DownFly/alns/",Files[i],".rda"))
return(aln)
}
stopCluster(cl)

MuscleToSequences = function(aln) { # define a function to process alingment data
aln = aln$seq # parse object 
sequences = list()
for(i in 1:nrow(aln)) {
x = aln[i,2] # grab second column with seq
sequences[[i]] = tolower(substring(x, seq(1,nchar(x),1), seq(1,nchar(x),1))) # split into ind characaters
}
names(sequences) = aln[,1] # grab the taxa names

return(sequences)
}

ProcessAln = function(X) {
SeqList = list()
for(i in 1:length(X)) {
sequences = X[[i]]
sequences = MuscleToSequences(sequences) # to list of sequences
# sequences = sample(sequences,ceiling(length(sequences)*0.7)) # randomly remove 70% for testing ONLY for testing 
sequences = sequences[!duplicated(sequences)] # remove duplicates sequences
names(sequences) = make.unique(names(sequences), sep = "") # make names unique
SeqList[[i]] = sequences
}

return(SeqList)
}


AlnList = ProcessAln(X)
names(AlnList) = gsub(".fasta","",Files)

save(AlnList, file = "C:/Users/Owner/Desktop/BodySize/DownFly/AlnList.rda")
}

if(FALSE) { # 7b. # Run muscle alignment NOT in parallel
Dir = "C:/Users/Owner/Desktop/BodySize/DownFly/Fastas/"
Files = list.files(Dir)
Done = list.files("C:/Users/Owner/Desktop/BodySize/DownFly/alns/")

Genes = c("myosin light chain light chain","tRNA-Val","histone H3","12S ribosomal RNA","NADH dehydogenase subunit I","tRNA-Leu","cytochrome oxidase II","elongation factor 1 alpha","internal transcribed spacer 2","5.8S ribosomal RNA","internal transcribed spacer 1","18S ribosomal RNA","28S ribosomal RNA","16S ribosomal RNA","cytochrome oxidase I")

Genes = Genes[!Genes %in% gsub(".rda|.fasta","",Done)] # don't do those already done
Files = Files[gsub(".fasta","",Files) %in% Genes] # remove files not in Genes

library(muscle)
for(i in 1:length(Files)) {
File = paste0(Dir,Files[i])
print(File)
aln = muscle(seqs = File)
save(aln,file = paste0("C:/Users/Owner/Desktop/BodySize/DownFly/alns/",Files[i],".rda"))
}

}

# Alignment Cleaning ------>

if(FALSE) { # 8. # Run Gblock alignment clean up 
# Be aware windows causes a false error sometimes. Download new version replace old. 

load(file = "C:/Users/Owner/Desktop/BodySize/DownFly/AlnList.rda")

library(ape)
library(phyloch)

gblocks = function (x, b1 = 0.5, b2 = b1, b3 = ncol(x), b4 = 2, b5 = "a", exec) { # personlized gblocks function
	if (inherits(x, "alignment")) 
        x <- as.DNAbin(x)
    if (inherits(x, "list")) 
        stop("cannot handle unaligned sequences")
    if (b1 < 0.5 | b1 > 1) 
        stop("b1 not in [0.5, 1]")
    if (b2 < b1 | b2 > 1) 
        stop("b2 not in [b1, 1]")
    if (b3 < 0 | b4 > ncol(x)) 
        stop("b3 not in [0, ", ncol(x), "]")
    if (b4 < 2 | b4 > ncol(x)) 
        stop("b4 not in [2, ", ncol(x), "]")
    b5 <- match.arg(b5, c("a", "h", "n"))
    rwd <- getwd()
    if (missing(exec)) 
        exec <- "/Applications/Gblocks_0.91b"
    setwd(exec)
    ntax <- nrow(x)
	b1 <- round(ntax * b1) + 1
    b2 <- round(ntax * b2) + 1
    cat("--- Executing Gblocks: ---")
    cat("\nMinimum number of sequences for a conserved position:", 
        b1)
    cat("\nMinimum number of sequences for a flank position:", 
        b2)
    cat("\nMaximum number of contiguous nonconserved positions:", 
        b3)
    cat("\nMinimum length of a block:", b4)
    cat("\nAllowed gap positions:", b5)
    ape::write.dna(x, file = paste0(exec,"R2GBLOCK.fas"),format = "fasta")
	system(paste("./Gblocks R2GBLOCK.fas -t=d", " -b1=", b1, 
        " -b2=", b2, " -b3=", b3, " -b4=", b4, " -b5=", b5, sep = ""), 
        show.output.on.console = TRUE)
    Sys.sleep(3)
    out = ape::read.FASTA(paste0(exec,"R2GBLOCK.fas-gb"))
    Sys.sleep(3)
    do.call(file.remove,list(list.files(pattern = "R2GBLOCK")))
    setwd(rwd)
	return(out)
}

print(names(AlnList)) # check genes

AlnList_gb = list()
for(i in 1:length(AlnList)) {
sequences = AlnList[[i]]
names(sequences) = gsub(" ","_",names(sequences)) # clean up names
x = do.call(rbind,sequences) # turn into a matrix to be used by gblocks
AlnList_gb[[i]] = as.character(gblocks(x, b1 = 0.5, b2 = 0.5, b3 = 8, b4 = 10, b5 = "a", exec = "C:/Users/Owner/Desktop/BodySize/DownFly/Gblocks_0.91b/")) # parameters b3 and b4 are the important ones here

}
names(AlnList_gb) = names(AlnList) # rename

save(AlnList_gb, file = "C:/Users/Owner/Desktop/BodySize/DownFly/AlnList_gb.rda")
}

if(FALSE) { # 9. # Run partition finder...

source(file = "C:/Users/Owner/Desktop/Scripts/R/TreeBuildPackage.r")
load(file = "C:/Users/Owner/Desktop/BodySize/DownFly/AlnList_gb.rda")

# Write cfg file
top=c(
"## ALIGNMENT FILE ##
alignment = test.phy;

## BRANCHLENGTHS: linked | unlinked ##
branchlengths = linked;

## MODELS OF EVOLUTION for PartitionFinder: all | raxml | mrbayes | beast | <list> ##
##              for PartitionFinderProtein: all_protein | <list> ##
models = all;

# MODEL SELECCTION: AIC | AICc | BIC #
model_selection = BIC;

## DATA BLOCKS: see manual for how to define ##
[data_blocks]\n")

Lengths = sapply(AlnList_gb,function(x) length(x[[1]]))

# Lengths
Start = c(1,unname(cumsum(Lengths) + 1)[2:length(Lengths)-1])
Finish = unname((cumsum(Lengths)))[2:length(Lengths)-1]
Finish = unname(c(Finish,tail(cumsum(Lengths),1)))

Names = gsub(" ","_",names(AlnList_gb))
Names = gsub(".","",Names) # replace strange characters
partitions = paste0(Names," = ",Start,"-",Finish,";",collapse="\n")

bottom = c("

## SCHEMES, search: all | greedy | rcluster | hcluster | user ##
[schemes]
search = greedy;

#user schemes go here if search=user. See manual for how to define.#")

x = paste0(top,partitions,bottom)
write(x, file = "C:/Users/Owner/Desktop/test.cfg",append = FALSE, sep = " ")


# Write data file
AlnList_gb = Filler(AlnList_gb) # fill genes
sequences = ConcatGenes(AlnList_gb) # combine genes

# write file formate for partition finder
nsp = length(sequences) # number of species
nns = length(sequences[[1]]) # number of nucleotides

D = data.frame(sp = names(sequences),seqs = sapply(sequences,function(x) toupper(paste(x,collapse=""))))
D = rbind(data.frame(sp=nsp,seqs=nns),D)

write.table(D,file="C:/Users/Owner/Desktop/BodySize/DownFly/PartitionFinder/downfly/downfly.phy",quote=FALSE,row.names=FALSE,col.names=FALSE,sep="    ")

# python PartitionFinder.py "C:/Users/Owner/Desktop/BodySize/DownFly/PartitionFinder/downfly/" --raxml # run this command using system
}

if(FALSE) { # 10. # Read best_scheme.txt and package into Partitions.rda
# we will use this to divide up our AlnList_gb.rda optimally
Lines = readLines("C:/Users/Owner/Desktop/BodySize/DownFly/PartitionFinder/downfly/analysis/best_scheme.txt")

# clean file
Lines = Lines[which(grepl("Subset | Best Model | Subset Partitions              | Subset Sites                   | Alignment                               ",Lines)):length(Lines)]
Lines = Lines[2:(which(Lines == "")[1]-1)] # first empty line
Lines = sapply(Lines,function(x) strsplit(x,"\\|") )
Lines = lapply(Lines,function(x) gsub(" ","",x)) # remove whitespace
Lines = unname(lapply(Lines,function(x) data.frame(model = x[2], product = x[3], partition = x[4], stringsAsFactors = FALSE))) # data.frames
Lines = do.call("rbind",Lines)
Lines$partnum = 1:nrow(Lines) # partition number

P = strsplit(Lines$partition,",") # Partitions split
L = list(); for(i in 1:length(P)) L[[i]] = data.frame(partnum = Lines$partnum[i],model = Lines$model[i],partition = P[[i]],product = Lines$product[i],stringsAsFactors = FALSE)
Lines = do.call("rbind",L)
sf = as.data.frame(do.call("rbind",strsplit(Lines$partition,"-"))) # start finish matrix
colnames(sf) = c("start","finish")

Partitions = cbind(Lines,sf) # got our partitions

# clean up
Partitions$start = as.numeric(as.character(Partitions$start))
Partitions$finish = as.numeric(as.character(Partitions$finish))

# save 
# save(Partitions,file="C:/Users/Owner/Desktop/BodySize/DownFly/Partitions.rda")
}

if(FALSE) { # 11. # divide up AlnList_gb.rda according to PartitionFinder to Alignments.rda"
source(file = "C:/Users/Owner/Desktop/Scripts/R/TreeBuildPackage.r")
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Partitions.rda")
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/AlnList_gb.rda")

# combine all genes into one alignment 
AlnList_gb = Filler(AlnList_gb) # fill genes
sequences = ConcatGenes(AlnList_gb) # combine genes

S = sequences 
P = split(Partitions,Partitions$partnum) # partitions list

grab = function(s,p) { # grab partition and combine it
GL = list()
for(i in 1:nrow(p)) {
GL[[i]] = s[p$start[i]:p$finish[i]] # grab list
}
grabbed = do.call("c",GL) # get grabbed partition
return(grabbed)
}

# do all the grabbing and combining 
out = list()
for(j in 1:length(S)) {
s = S[[j]] # focal sequence
sp = names(S)[j] # focal sp

part_list = list()
for(i in 1:length(P)) {
p = P[[i]]
part_list[[i]] = grab(s,p) 
}

pl = list(part_list) # should be of length 11
names(pl) = sp
out[[j]] = pl
}

Alignments = out # get alignments for Beast
# names(Alignments) = names(S) # rename

save(Alignments,file="C:/Users/Owner/Desktop/BodySize/DownFly/Alignments.rda") # 
}

# Run in BEAST --------> 

if(FALSE) { # Zygoptera
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Alignments.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Species.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Partitions.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Fossils.rda") 

# Species
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies)

Outgroup = Species[Species$SubOrder == "Anisoptera",] # get dameselflies for outgroup
Species = Species[Species$SubOrder == "Zygoptera",]

# Alignments
A = Alignments 
A = sapply(A,function(x) x) # go up one level 
taxa_with = names(A) # get taxa with genes

# Get taxa with genes
outgroup_with = taxa_with[taxa_with %in% Outgroup$GenusSpecies] 
taxa_with = taxa_with[taxa_with %in% Species$GenusSpecies] 

# Subset taxa into those only with genes
Species = Species[Species$GenusSpecies %in% taxa_with,] # only those species with genes
Outgroup = Outgroup[Outgroup$GenusSpecies %in% outgroup_with,] # only those species with genes
Outgroup = Outgroup[sample(1:nrow(Outgroup),50),] # randomly draw from outgroup

# Taxon sets # T
T_Genus = lapply(split(Species,Species$Genus),function(x) x$GenusSpecies) # split into groups
T_Family = lapply(split(Species,Species$Family),function(x) x$GenusSpecies) # split into groups

T = c(T_Genus,T_Family)

T$outgroup = Outgroup$GenusSpecies
T$ingroup = Species$GenusSpecies
T$taxa = c(Species$GenusSpecies,Outgroup$GenusSpecies)

T = T[!sapply(T,length) <= 1] # remove those groups with just one

# Alignments # A
A = A[names(A) %in% c(Species$GenusSpecies,Outgroup$GenusSpecies)] 
AL = list(); for(i in 1:length(A[[1]])) AL[[i]] = lapply(A,function(x) x[[i]]); A = AL # unlist, reorder

a = paste0("alignment",1:length(A))
x = names(T) # taxonset names

# mcmc
chainLength = "10 000 000"
operatorAnalysis = "Dragonfly.ops.txt"

# Fossils # F

GetFossils = function(Genera = "Sapho", Fossils) { # Get Genera with fossils

FossilsExtractor = function(Genus,Fossils) { # define fossils function 
# max(Fossils[Fossils$Genus == Genus,]$ma_mid)
max(Fossils[Fossils$Genus == Genus,]$ma_max)
}

FL = list() # fossil list
for(i in 1:length(Genera)){
FL[[i]] = c(x = Genera[i],stdev = 2,Mean=FossilsExtractor(Genus = Genera[i],Fossils))
}

return(FL)
}

TaxaWithFossils = unique(Species$Genus[Species$Genus %in% Fossils$Genus]) 
F = GetFossils(Genera = TaxaWithFossils, Fossils) # should be 250
F = F[!sapply(F,function(x) x["x"]) %in% names(table(Species$Genus)[table(Species$Genus) <= 1])] # remove those groups with only one species
R = c(x = "outgroup",stdev = "3", Mean = "132.9")

# Ware Fossils "How to date a Dragonfly 2016"
# Odonata 237 - root
# Zygotera 132.9
# Epiophlebia (Anisozygoptera) + Anisoptera 199.3
# Anisoptera 168
# Aeshnidae 139.8
# Gomphidae 150
# Libellulidae 29.2 
# Macromiidae 15.5
# Corduliidae 12.7 # remove Cordulia

# Ware "Going with the flow" inferred dates 


StartTree = function(D,Outgroup) { # start tree function
library(ape)
D$GenusSpecies = gsub(" ","_",D$GenusSpecies)
D = as.data.frame(unclass(D))
D = rbind(D,Outgroup)
tree = as.phylo(~SubOrder/Family/Genus/GenusSpecies, data=D)
tree = root(tree,Outgroup$GenusSpecies,resolve.root=TRUE)
return(tree)
}

StartTree = StartTree(D = Species,Outgroup)
StartTree$edge.length = rep(10, nrow(StartTree$edge)) # add tip lengths
StartTree = write.tree(StartTree) # make newick tree
File = "Zygoptera"
}

if(FALSE) { # Anisoptera
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Alignments.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Species.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Partitions.rda") 
load(file="C:/Users/Owner/Desktop/BodySize/DownFly/Fossils.rda") 

# Species
Species$GenusSpecies = gsub(" ","_",Species$GenusSpecies)

Outgroup = Species[Species$SubOrder == "Zygoptera",] # get dameselflies for outgroup
Species = Species[Species$SubOrder == "Anisoptera",]

# Alignments
A = Alignments 
A = sapply(A,function(x) x) # go up one level 
taxa_with = names(A) # get taxa with genes

# Get taxa with genes
outgroup_with = taxa_with[taxa_with %in% Outgroup$GenusSpecies] 
taxa_with = taxa_with[taxa_with %in% Species$GenusSpecies] 

# Subset taxa into those only with genes
Species = Species[Species$GenusSpecies %in% taxa_with,] # only those species with genes
Outgroup = Outgroup[Outgroup$GenusSpecies %in% outgroup_with,] # only those species with genes
Outgroup = Outgroup[sample(1:nrow(Outgroup),50),] # randomly draw from outgroup

# Taxon sets # T
T_Genus = lapply(split(Species,Species$Genus),function(x) x$GenusSpecies) # split into groups
T_Family = lapply(split(Species,Species$Family),function(x) x$GenusSpecies) # split into groups

T = c(T_Genus,T_Family)

T$outgroup = Outgroup$GenusSpecies
T$ingroup = Species$GenusSpecies
T$taxa = c(Species$GenusSpecies,Outgroup$GenusSpecies)

T = T[!sapply(T,length) <= 1] # remove those groups with just one

# Alignments # A
A = A[names(A) %in% c(Species$GenusSpecies,Outgroup$GenusSpecies)] 
AL = list(); for(i in 1:length(A[[1]])) AL[[i]] = lapply(A,function(x) x[[i]]); A = AL # unlist, reorder

a = paste0("alignment",1:length(A))
x = names(T) # taxonset names

# mcmc
chainLength = "10 000 000"
operatorAnalysis = "Anisoptera.ops.txt"

# Fossils # F

GetFossils = function(Genera = "Sapho", Fossils) { # Get Genera with fossils

FossilsExtractor = function(Genus,Fossils) { # define fossils function 
# max(Fossils[Fossils$Genus == Genus,]$ma_mid)
max(Fossils[Fossils$Genus == Genus,]$ma_max)
}

FL = list() # fossil list
for(i in 1:length(Genera)){
FL[[i]] = c(x = Genera[i],stdev = 2,Mean=FossilsExtractor(Genus = Genera[i],Fossils))
}

return(FL)
}

TaxaWithFossils = unique(Species$Genus[Species$Genus %in% Fossils$Genus]) 
F = GetFossils(Genera = TaxaWithFossils, Fossils) # should be 250
F = F[!sapply(F,function(x) x["x"]) %in% names(table(Species$Genus)[table(Species$Genus) <= 1])] # remove those groups with only one species
R = c(x = "outgroup",stdev = "3", Mean = "230")

# add and remove fossils base on "How to date a Dragonfly" Ware 2016

# Ware Fossils "How to date a Dragonfly 2016"
# Odonata 237 - root in R
# Zygotera 132.9 x
# Epiophlebia (Anisozygoptera) + Anisoptera 199.3 x
# Anisoptera 168 # ingroup x
# Aeshnidae 139.8 yes
# Gomphidae 150 yes
# Libellulidae 29.2 x 
# Macromiidae 15.5 x
# Corduliidae 12.7 yes # Cordulia has suspect old fossil

F = F[!sapply(F,function(x) x[1] == "Cordulia")] # remove cordulia because of suspect fossil

F_add = list(
c(x="Aeshnidae",stdev="2",Mean="139.8"),
c(x="ingroup",stdev="2",Mean="168"),
c(x="Gomphidae",stdev="2",Mean="150"),
c(x="Corduliidae",stdev="2",Mean="12.7")
)

F = c(F,F_add) # add missing dates

StartTree = function(D,Outgroup) { # start tree function
library(ape)
D$GenusSpecies = gsub(" ","_",D$GenusSpecies)
D = as.data.frame(unclass(D))
D = rbind(D,Outgroup)
tree = as.phylo(~SubOrder/Family/Genus/GenusSpecies, data=D)
tree = root(tree,Outgroup$GenusSpecies,resolve.root=TRUE)
return(tree)
}

StartTree = StartTree(D = Species,Outgroup)
StartTree$edge.length = rep(10, nrow(StartTree$edge)) # add tip lengths
StartTree = write.tree(StartTree) # make newick tree
File = "Anisoptera"

# run in command
# java -Xms64m -Xmx4000m -jar lib/beast.jar 
}

if(FALSE) { # BEAST XML

# key variables 
# T: taxa sets of genera, families, ingroup ect 
# x: names from T

# A: alignment list
# a: names from alignment list, alignment1, alignment2 ...

# F: Fossils list
# R: root list


# beast
beast = '<beast>\n'

TaxonSets = function(T) { # taxa

TaxaSet = function(T) {
xtaxa = T$taxa

q = "\"" # quote mark
Open = paste0('<taxa id="taxa">\n',collapse="")
Meat = paste0("<taxon id=",q,xtaxa,q,"/>",collapse="\n")
Close = "\n</taxa>\n"

out = paste0(Open,Meat,Close,collapse="\n")
return(out)
}

out_taxa = TaxaSet(T)

T$taxa = NULL # remove taxa

TaxonSet = function(x,set) {

q = "\"" # quote mark
Open = paste0("<taxa id=",q,set,q,">\n",collapse="")
Meat = paste0("<taxon idref=",q,x,q,"/>",collapse="\n")
Close = "\n</taxa>\n"

out = paste0(Open,Meat,Close,collapse="\n")
return(out)
}

out = list()
for(i in 1:length(T)) out[[i]] = TaxonSet(T[[i]],names(T)[i]) 
out = paste0(out,collapse="\n")

out = paste0(out_taxa,out,collapse="\n")
return(out)
}

taxa = TaxonSets(T) 

AlignmentSets = function(A) { # alignment

AlignmentSet = function(a,set) {
a = lapply(a,function(x) toupper(paste(x,collapse="")))

q = "\"" # quote mark

Open = paste0('<alignment id=',q,set,q,' dataType="nucleotide">\n',collapse = "")
Meat = paste0(
'<sequence>\n',
'<taxon idref=',q,names(a),q,'/>\n',
unlist(a),"\n",
'</sequence>\n',
collapse="")
Close = '</alignment>\n'

out = paste0(Open,Meat,Close,collapse="\n")
return(out)
}

out = list()
for(i in 1:length(A)) {
set = paste0("alignment",i)
a = A[[i]]
out[[i]] = AlignmentSet(a,set)
}

out = paste0(out,collapse="\n")
return(out)
}

alignment = AlignmentSets(A)

PatternSets = function(a) { # patterns

q = "\"" # quote mark
out = paste0(
'<patterns id=',q,a,'.patterns',q,' from="1" strip="false">\n',
'<alignment idref=',q,a,q,'/>\n',
'</patterns>\n',
collapse="\n")

return(out)
}

patterns = PatternSets(a)

YuleModel = function() { # yuleModel

out = paste0(
'<yuleModel id="yule" units="years">\n',
'<birthRate>\n',
'<parameter id="yule.birthRate" value="3.6" lower="0.0"/>\n',
'</birthRate>\n',
'</yuleModel>\n',
collapse="\n")

return(out)
}

yuleModel = YuleModel()

ConstantSize = function() { # constantSize

out = paste0(
'<constantSize id="initialDemo" units="years">\n',
'<populationSize>\n',
'<parameter id="initialDemo.popSize" value="100.0"/>\n',
'</populationSize>\n',
'</constantSize>\n',
collapse="\n")

return(out)
}

constantSize = ConstantSize()

CoalescentSimulator = function(x) { # coalescentSimulator
x = x[!x %in% c("ingroup","outgroup","taxa")] # remove outgroup and taxa

# filter outgroup and taxa

q = "\"" # quote mark
Open = paste0(
'<coalescentSimulator id="startingTree">\n',
'<coalescentSimulator>\n',
collapse="\n")

Meat = paste0(
'<coalescentSimulator>\n',
'<taxa idref=',q,x,q,'/>\n',
'<constantSize idref="initialDemo"/>\n',
'</coalescentSimulator>\n',
collapse="\n")

Close = paste0(
'<taxa idref="ingroup"/>\n',
'<constantSize idref="initialDemo"/>\n',
'</coalescentSimulator>\n',
'<coalescentSimulator>\n',
'<taxa idref="outgroup"/>\n',
'<constantSize idref="initialDemo"/>\n',
'</coalescentSimulator>\n',
'<taxa idref="taxa"/>\n',
'<constantSize idref="initialDemo"/>\n',
'</coalescentSimulator>\n',
collapse="\n")

out = paste0(Open,Meat,Close,collapse="\n")

return(out)
}

# coalescentSimulator = CoalescentSimulator(x)


TreeModel = function(StartTree) { # treeModel
out = paste0(
'<newick id="startingTree">',StartTree,'</newick>\n',
'<treeModel id="treeModel">\n',
'<newick idref="startingTree"/>\n',
'<rootHeight>\n',
'<parameter id="treeModel.rootHeight"/>\n',
'</rootHeight>\n',
'<nodeHeights internalNodes="true">\n',
'<parameter id="treeModel.internalNodeHeights"/>\n',
'</nodeHeights>\n',
'<nodeHeights internalNodes="true" rootNode="true">\n',
'<parameter id="treeModel.allInternalNodeHeights"/>\n',
'</nodeHeights>\n',
'</treeModel>\n',
collapse="\n")

return(out)
}

treeModel = TreeModel(StartTree)

TmrcaStatistic = function(x) { # tmrcaStatistic
x = x[!x == "taxa"] # remove taxa part
# x = x[!x == "outgroup"] # remove outgroup

q = "\"" # quote mark
out = paste0(
'<tmrcaStatistic id=',q,'tmrca(',x,')',q,'>\n',
'<mrca>\n',
'<taxa idref=',q,x,q,'/>\n',
'</mrca>\n',
'<treeModel idref="treeModel"/>\n',
'</tmrcaStatistic>\n',
'<monophylyStatistic id=',q,'monophyly(',x,')',q,'>\n',
'<mrca>\n',
'<taxa idref=',q,x,q,'/>\n',
'</mrca>\n',
'<treeModel idref="treeModel"/>\n',
'</monophylyStatistic>\n',
collapse="\n")

return(out)
}

tmrcaStatistic = TmrcaStatistic(x)

SpeciationLikelihood = function() { # speciationLikelihood

out = paste0('<speciationLikelihood id="speciation">\n',
'<model>\n',
'<yuleModel idref="yule"/>\n',
'</model>\n',
'<speciesTree>\n',
'<treeModel idref="treeModel"/>\n',
'</speciesTree>\n',
'</speciationLikelihood>\n',
collapse="\n")

return(out)
}

speciationLikelihood = SpeciationLikelihood()

DiscretizedBranchRates = function() { # discretizedBranchRates

out = paste0(
'<discretizedBranchRates id="branchRates">\n',
'<treeModel idref="treeModel"/>\n',
'<distribution>\n',
'<logNormalDistributionModel meanInRealSpace="true">\n',
'<mean>\n',
'<parameter id="ucld.mean" value="1.0" lower="0.0"/>\n',
'</mean>\n',
'<stdev>\n',
'<parameter id="ucld.stdev" value="0.3333333333333333" lower="0.0"/>\n',
'</stdev>\n',
'</logNormalDistributionModel>\n',
'</distribution>\n',
'<rateCategories>\n',
'<parameter id="branchRates.categories"/>\n',
'</rateCategories>\n',
'</discretizedBranchRates>\n',
collapse="\n")

return(out)
}

discretizedBranchRates = DiscretizedBranchRates()

RateStatistic = function() { # rateStatistic
out = paste0(
'<rateStatistic id="meanRate" name="meanRate" mode="mean" internal="true" external="true">\n',
'<treeModel idref="treeModel"/>\n',
'<discretizedBranchRates idref="branchRates"/>\n',
'</rateStatistic>\n',
'<rateStatistic id="coefficientOfVariation" name="coefficientOfVariation" mode="coefficientOfVariation" internal="true" external="true">\n',
'<treeModel idref="treeModel"/>\n',
'<discretizedBranchRates idref="branchRates"/>\n',
'</rateStatistic>\n',
'<rateCovarianceStatistic id="covariance" name="covariance">\n',
'<treeModel idref="treeModel"/>\n',
'<discretizedBranchRates idref="branchRates"/>\n',
'</rateCovarianceStatistic>\n',
collapse="\n")

return(out)
}

rateStatistic = RateStatistic()

GtrModel = function(a) { # gtrModel

q = "\"" # quote mark
out = paste0('<gtrModel id=',q,a,'.gtr">\n',
'<frequencies>\n',
'<frequencyModel dataType="nucleotide">\n',
'<frequencies>\n',
'<parameter id=',q,a,'.frequencies" value="0.25 0.25 0.25 0.25"/>\n',
'</frequencies>\n',
'</frequencyModel>\n',
'</frequencies>\n',
'<rateAC>\n',
'<parameter id=',q,a,'.ac" value="1.0" lower="0.0"/>\n',
'</rateAC>\n',
'<rateAG>\n',
'<parameter id=',q,a,'.ag" value="1.0" lower="0.0"/>\n',
'</rateAG>\n',
'<rateAT>\n',
'<parameter id=',q,a,'.at" value="1.0" lower="0.0"/>\n',
'</rateAT>\n',
'<rateCG>\n',
'<parameter id=',q,a,'.cg" value="1.0" lower="0.0"/>\n',
'</rateCG>\n',
'<rateGT>\n',
'<parameter id=',q,a,'.gt" value="1.0" lower="0.0"/>\n',
'</rateGT>\n',
'</gtrModel>\n',
'<siteModel id=',q,a,'.siteModel">\n',
'<substitutionModel>\n',
'<gtrModel idref=',q,a,'.gtr"/>\n',
'</substitutionModel>\n',
'<gammaShape gammaCategories="4">\n',
'<parameter id=',q,a,'.alpha" value="0.5" lower="0.0"/>\n',
'</gammaShape>\n',
'<proportionInvariant>\n',
'<parameter id=',q,a,'.pInv" value="0.5" lower="0.0" upper="1.0"/>\n',
'</proportionInvariant>\n',
'</siteModel>\n',
collapse="\n")

return(out)
}

TreeLikelihood1 = function(a) { # treeLikelihood

q = "\"" # quote mark
out = paste0(
'<treeLikelihood id=',q,a,'.treeLikelihood" useAmbiguities="false">\n',
'<patterns idref=',q,a,'.patterns"/>\n',
'<treeModel idref="treeModel"/>\n',
'<siteModel idref=',q,a,'.siteModel"/>\n',
'<discretizedBranchRates idref="branchRates"/>\n',
'</treeLikelihood>\n',
collapse="\n")

return(out)
}

treeLikelihood1 = TreeLikelihood1(a)

gtrModel = GtrModel(a)

ScaleOperator = function(a) { # scaleOperator

q = "\"" # quote mark
Open = '<operators id="operators" optimizationSchedule="default">\n'

Meat = paste0(
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.ac"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.ag"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.at"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.cg"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.gt"/>\n',
'</scaleOperator>\n',
'<deltaExchange delta="0.01" weight="0.1">\n',
'<parameter idref=',q,a,'.frequencies"/>\n',
'</deltaExchange>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.alpha"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="0.1">\n',
'<parameter idref=',q,a,'.pInv"/>\n',
'</scaleOperator>\n',
collapse="\n")

Close = paste0(
'<scaleOperator scaleFactor="0.75" weight="3">\n',
'<parameter idref="ucld.mean"/>\n',
'</scaleOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="3">\n',
'<parameter idref="ucld.stdev"/>\n',
'</scaleOperator>\n',
'<subtreeSlide size="1.1" gaussian="true" weight="15">\n',
'<treeModel idref="treeModel"/>\n',
'</subtreeSlide>\n',
'<narrowExchange weight="15">\n',
'<treeModel idref="treeModel"/>\n',
'</narrowExchange>\n',
'<wideExchange weight="3">\n',
'<treeModel idref="treeModel"/>\n',
'</wideExchange>\n',
'<wilsonBalding weight="3">\n',
'<treeModel idref="treeModel"/>\n',
'</wilsonBalding>\n',
'<scaleOperator scaleFactor="0.75" weight="3">\n',
'<parameter idref="treeModel.rootHeight"/>\n',
'</scaleOperator>\n',
'<uniformOperator weight="30">\n',
'<parameter idref="treeModel.internalNodeHeights"/>\n',
'</uniformOperator>\n',
'<scaleOperator scaleFactor="0.75" weight="3">\n',
'<parameter idref="yule.birthRate"/>\n',
'</scaleOperator>\n',
'<upDownOperator scaleFactor="0.75" weight="3">\n',
'<up>\n',
'<parameter idref="ucld.mean"/>\n',
'</up>\n',
'<down>\n',
'<parameter idref="treeModel.allInternalNodeHeights"/>\n',
'</down>\n',
'</upDownOperator>\n',
'<swapOperator size="1" weight="10" autoOptimize="false">\n',
'<parameter idref="branchRates.categories"/>\n',
'</swapOperator>\n',
'<uniformIntegerOperator weight="10">\n',
'<parameter idref="branchRates.categories"/>\n',
'</uniformIntegerOperator>\n',
'</operators>\n',
collapse="\n")

out = paste0(Open,Meat,Close,collapse="\n")
return(out)
}

scaleOperator = ScaleOperator(a)

MCMC = function(chainLength,operatorAnalysis) { # mcmc
chainLength = gsub(" ","",chainLength)
q = "\"" # quote mark
out = paste0(
'<mcmc id="mcmc" chainLength=',q,chainLength,q,' autoOptimize="true" operatorAnalysis=',q,operatorAnalysis,q,'>\n',
'<posterior id="posterior">\n',
'<prior id="prior">\n',
collapse="\n")
return(out)
}

mcmc = MCMC(chainLength,operatorAnalysis)

BooleanLikelihood = function(x) { # booleanLikelihood
x = x[!x == "taxa"] # remove taxa
q = "\"" # quote mark
Open = '<booleanLikelihood>\n'
Meat = paste0('<monophylyStatistic idref="monophyly(',x,')"/>\n',collapse="\n")
Close = '</booleanLikelihood>\n'

out = paste0(Open,Meat,Close,collapse="\n")

return(out)
}

booleanLikelihood = BooleanLikelihood(x)

LogNormalPrior = function(R,F) { # logNormalPrior.............. # Fossils Calibrations

q = "\"" # quote mark
lognorm = list()
for(i in 1:length(F)) {
f = F[[i]]
lognorm[[i]] = paste0(
'<logNormalPrior mean=',q,f["Mean"],q,' stdev=',q,f["stdev"],q,' offset="0.0" meanInRealSpace="true">\n',
'<statistic idref="tmrca(',f["x"],')"/>\n',
'</logNormalPrior>\n',
collapse="\n")
}

out1 = paste0(do.call("c",lognorm),collapse="\n")

out2 = paste0(
'<normalPrior mean=',q,R["Mean"],q,' stdev=',q,R["stdev"],q,'>\n',
'<statistic idref="tmrca(',R["x"],')"/>\n',
'</normalPrior>\n',
collapse="\n")

out = paste0(out1,out2,collapse="\n")

return(out)
}

logNormalPrior = LogNormalPrior(R,F)

GammaPrior = function(a) { # gammaPrior # all boilerplate priors

q = "\"" # quote mark
Meat = paste0(
'<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n',
'<parameter idref=',q,a,'.ac"/>\n',
'</gammaPrior>\n',
'<gammaPrior shape="0.05" scale="20.0" offset="0.0">\n',
'<parameter idref=',q,a,'.ag"/>\n',
'</gammaPrior>\n',
'<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n',
'<parameter idref=',q,a,'.at"/>\n',
'</gammaPrior>\n',
'<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n',
'<parameter idref=',q,a,'.cg"/>\n',
'</gammaPrior>\n',
'<gammaPrior shape="0.05" scale="10.0" offset="0.0">\n',
'<parameter idref=',q,a,'.gt"/>\n',
'</gammaPrior>\n',
'<uniformPrior lower="0.0" upper="1.0">\n',
'<parameter idref=',q,a,'.frequencies"/>\n',
'</uniformPrior>\n',
'<exponentialPrior mean="0.5" offset="0.0">\n',
'<parameter idref=',q,a,'.alpha"/>\n',
'</exponentialPrior>\n',
'<uniformPrior lower="0.0" upper="1.0">\n',
'<parameter idref=',q,a,'.pInv"/>\n',
'</uniformPrior>\n',
collapse="\n")

Close = 
paste0(
'<exponentialPrior mean="0.3333333333333333" offset="0.0">\n',
'<parameter idref="ucld.stdev"/>\n',
'</exponentialPrior>\n',
'<uniformPrior lower="0.0" upper="1.0E100">\n',
'<parameter idref="yule.birthRate"/>\n',
'</uniformPrior>\n',
'<speciationLikelihood idref="speciation"/>\n',
'</prior>\n',
collapse="\n")

out = paste0(Meat,Close,collapse="\n")

return(out)
}

gammaPrior = GammaPrior(a)

TreeLikelihood2 = function(a) { # treeLikelihood

q = "\"" # quote mark
Open = '<likelihood id="likelihood">\n'
Meat = paste0(
'<treeLikelihood idref=',q,a,'.treeLikelihood"/>\n',
collapse="\n")
Close = paste0(
'</likelihood>\n',
'</posterior>\n',
'<operators idref="operators"/>\n',
collapse="\n")

out = paste0(Open,Meat,Close,collapse="\n")

return(out)
}

treeLikelihood2 = TreeLikelihood2(a)

ScreenLog = function() { # screenLog # boilerplate screelog xml

out = paste0(
'<log id="screenLog" logEvery="1000">\n',
'<column label="Posterior" dp="4" width="12">\n',
'<posterior idref="posterior"/>\n',
'</column>\n',
'<column label="Prior" dp="4" width="12">\n',
'<prior idref="prior"/>\n',
'</column>\n',
'<column label="Likelihood" dp="4" width="12">\n',
'<likelihood idref="likelihood"/>\n',
'</column>\n',
'<column label="rootHeight" sf="6" width="12">\n',
'<parameter idref="treeModel.rootHeight"/>\n',
'</column>\n',
'<column label="ucld.mean" sf="6" width="12">\n',
'<parameter idref="ucld.mean"/>\n',
'</column>\n',
'</log>\n',
collapse="\n")

return(out)
}

screenLog = ScreenLog()

FileLog = function(a,x) { # fileLog
x = x[!x == "taxa"] # remove taxa
q = "\"" # quote mark
Open = paste0(
'<log id="fileLog" logEvery="1000" fileName=',q,File,'.log.txt" overwrite="false">\n',
'<posterior idref="posterior"/>\n',
'<prior idref="prior"/>\n',
'<likelihood idref="likelihood"/>\n',
'<parameter idref="treeModel.rootHeight"/>\n',
'<parameter idref="yule.birthRate"/>\n',
'<parameter idref="ucld.mean"/>\n',
'<parameter idref="ucld.stdev"/>\n',
'<rateStatistic idref="meanRate"/>\n',
'<rateStatistic idref="coefficientOfVariation"/>\n',
'<rateCovarianceStatistic idref="covariance"/>\n',
'<speciationLikelihood idref="speciation"/>\n',
collapse="\n")

Meat1 = paste0(
'<tmrcaStatistic idref="tmrca(',x,')"/>\n',
collapse="\n")

Meat2 = paste0(
'<parameter idref=',q,a,'.ac"/>\n',
'<parameter idref=',q,a,'.ag"/>\n',
'<parameter idref=',q,a,'.at"/>\n',
'<parameter idref=',q,a,'.cg"/>\n',
'<parameter idref=',q,a,'.gt"/>\n',
'<parameter idref=',q,a,'.frequencies"/>\n',
'<parameter idref=',q,a,'.alpha"/>\n',
'<parameter idref=',q,a,'.pInv"/>\n',
collapse="\n")

Meat3 = paste0(
'<treeLikelihood idref=',q,a,'.treeLikelihood"/>\n',
collapse="\n")

Close = '</log>\n'

out = paste0(Open,Meat1,Meat2,Meat3,Close,collapse="\n")

return(out)
}

fileLog = FileLog(a,x)

LogTree = function() { # logTree # boilerplate logging and trees.txt output and closing tags

q = "\"" # quote mark

out = paste0(
'<logTree id="treeFileLog" logEvery="1000" nexusFormat="true" fileName=',q,File,'.trees.txt" sortTranslationTable="true">\n',
'<treeModel idref="treeModel"/>\n',
'<trait name="rate" tag="rate">\n',
'<discretizedBranchRates idref="branchRates"/>\n',
'</trait>\n',
'<posterior idref="posterior"/>\n',
'</logTree>\n',
'</mcmc>\n',
'<report>\n',
'<property name="timer">\n',
'<mcmc idref="mcmc"/>\n',
'</property>\n',
'</report>\n',
'</beast>\n',
collapse="\n")

return(out)
}

logTree = LogTree()
# coalescentSimulator,

out = paste0(
beast,
taxa,
alignment,
patterns,
yuleModel,
constantSize,
treeModel,
tmrcaStatistic,
speciationLikelihood,
discretizedBranchRates,
rateStatistic,
gtrModel,
treeLikelihood1,
scaleOperator,
mcmc,
booleanLikelihood,
logNormalPrior,
gammaPrior,
treeLikelihood2,
screenLog,
fileLog,
logTree,
collapse="\n")

write(out, file = paste0("C:/Users/Owner/Desktop/BodySize/DownFly/",File,".xml"),append = FALSE, sep = " ")


}

if(FALSE) { # read in tree and save as tree.rda
# Process file produced by TreeAnnotator
library(ips)
tree = read.beast("C:/Users/Owner/Desktop/BodySize/DownFly/ZygopteraOutputTree.txt")
save(tree,file="C:/Users/Owner/Desktop/BodySize/DownFly/Zygoptera_tree.rda")

tree = read.beast("C:/Users/Owner/Desktop/BodySize/DownFly/AnisopteraOutputTree.txt")
save(tree,file="C:/Users/Owner/Desktop/BodySize/DownFly/Anisoptera_tree.rda")

}



