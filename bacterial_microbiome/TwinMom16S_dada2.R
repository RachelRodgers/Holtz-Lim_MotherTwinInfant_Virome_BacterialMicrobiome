#!/usr/bin/env Rscript

# TwinMom16S_dada2.R

#----- Load Libraries & Data -----#
library("ShortRead")
library("dada2")
library("ggplot2")
library("msa")
library("phangorn")
library("phyloseq")
library("tidyverse")

names(filteredFiles) <- sampleNames

#----- Dereplication -----#
# Create a list that will hold dereplication objects for each sample
singles <- vector("list", length(sampleNames))
names(singles) <- sampleNames

# Populate the list
for(sample in sampleNames) {
  derepF <- derepFastq(filteredFiles[[sample]])
  singles[[sample]] <- dada(derepF, err = errF, multithread = TRUE)
}

rm(derepF)

#----- Construct Sequence Table -----#
# Construct sequence table and remove chimeras
sequenceTable <- makeSequenceTable(singles)
sequenceTableNoChimeras <- removeBimeraDenovo(sequenceTable,  multithread = TRUE)
saveRDS(sequenceTableNoChimeras, "sequenceTableNoChimeras.RDS")

#----- Assign Taxonomy -----#
taxaRankNamesFull <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxaRankNamesTrunc <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

# GreenGenes
taxaGG <- assignTaxonomy(sequenceTableNoChimeras, ggDB, multithread = TRUE)
colnames(taxaGG) <- taxaRankNamesFull
# remove prefixes from taxonomy
taxaGGModified <- data.frame(taxaGG)
taxaGGModified <- rownames_to_column(taxaGGModified)
taxaGGModified <- map_dfc(taxaGGModified, 
                          function(x) {str_remove(x, pattern = "^[:lower:]__")})
taxaGGModified <- as.matrix(column_to_rownames(taxaGGModified))

# RDP
taxaRDP <- assignTaxonomy(sequenceTableNoChimeras, rdpDB, multithread = TRUE)
colnames(taxaRDP) <- taxaRankNamesTrunc
taxaRDPPlus <- addSpecies(taxaRDP, rdpDBSpecies)

# Silva
taxaSilva <- assignTaxonomy(sequenceTableNoChimeras, silvaDB, multithread = TRUE)
colnames(taxaSilva) <- taxaRankNamesTrunc
taxaSilvaPlus <- addSpecies(taxaSilva, silvaDBSpecies)

# HitDB
taxaHitDB <- assignTaxonomy(sequenceTableNoChimeras, hitDB, multithread = TRUE)

save.image("TwinMom16S_dada2_Pathogen_Server.RData") # read back w/ load() or attach()

#----- Construct Phylogenetic Tree -----#

# Get the sequences from the sequence table
seqs <- getSequences(sequenceTableNoChimeras)
names(seqs) <- seqs
# Multiple sequence alignment 
# (msa() is a long-running function that may take hours)
mult <- msa(seqs, method = "ClustalW", type = "dna", order = "input")
save.image("TwinMom16S_dada2_Pathogen_Server.RData")

# Convert MSA to phyDAT format
phangAlign <- as.phyDat(mult, type = "DNA",
                        names = getSequence(seqtab.nochim))
# Compute pairwise distances on phangAlign
dm <- dist.ml(phangAlign)
treeNJ <- NJ(dm)
# Compute likelihood of tree
fit <- pml(tree = treeNJ, data = phangAlign)
fitGTR <- update(fit, k = 4, inv = 0.2)
# (optim.pml() is a long-running function)
fitGTR <- optim.pml(fitGTR, model = "GTR", optInv = TRUE, optGamma = TRUE,
                    rearrangement = "stochastic", 
                    control = pml.control(trace = 0))

save.image("TwinMom16S_dada2_Pathogen_Server.RData") # read back w/ load() or attach()

#----- Build Phyloseq Objects -----#

# Greengenes
ps0.gg <- phyloseq(otu_table(sequenceTableNoChimeras, taxa_are_rows = FALSE),
                   tax_table(taxaGGModified), phy_tree(fitGTR$tree))
saveRDS(ps0.gg, "../data/physeqObjects/ps0.gg_single.RDS")

# RDP
ps0.rdp <- phyloseq(otu_table(sequenceTableNoChimeras, taxa_are_rows = FALSE),
                    tax_table(taxaRDPPlus), phy_tree(fitGTR$tree))
saveRDS(ps0.rdp, "../data/physeqObjects/ps0.rdp_single.RDS")

# silva
ps0.silva <- phyloseq(otu_table(sequenceTableNoChimeras, taxa_are_rows = FALSE),
                      tax_table(taxaSilvaPlus), phy_tree(fitGTR$tree))
saveRDS(ps0.silva, "../data/physeqObjects/ps0.silva_single.RDS")

# hitDB
ps0.hitdb <- phyloseq(otu_table(sequenceTableNoChimeras, taxa_are_rows = FALSE),
                      tax_table(taxaHitDB), phy_tree(fitGTR$tree))
saveRDS(ps0.hitdb, "../data/physeqObjects/ps0.hitdb_single.RDS")

#----- Save Data -----#

save.image("TwinMom16S_dada2_Pathogen_Server.RData")
writeLines(capture.output(sessionInfo()), 
           "TwinMom16S_dada2_Pathogen_Server_session_info.txt")
Sys.Date()
getwd()
sessionInfo()

