### Scripts used in Microbiome analysis

update.packages()
getwd()
setwd('/Users/benjaminwu/Dropbox/R.shared/R.1.Spain')

### Load Phyloseq
library(phyloseq)
library(ade4)
library(vegan)
library(biomformat)
library(devtools)
library(readr)
library(readtext)
library(tidyverse)
library(qiime2R)

install.packages('devtools')
install.packages('tidyverse')
install.packages('readr')
install.packages('readtext')
devtools::install_github("jbisanz/qiime2R")

### QIIME convert to JSON ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
biom convert -i /Users/benjaminwu/Dropbox/R.shared/R.1.Spain/Merged_otu_table.Nasal.Baseline.biom -o /Users/benjaminwu/Dropbox/R.shared/R.1.Spain/Merged_otu_table.Nasal.Baseline.txt --header-key "taxonomy" --to-tsv

biom convert -i /Users/benjaminwu/Dropbox/R.shared/R.1.Spain/Merged_otu_table.Nasal.Baseline.txt -o /Users/benjaminwu/Dropbox/R.shared/R.1.Spain/Merged_otu_table.Nasal.Baseline.json.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy
### QIIME convert to JSON ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Load the files needed, if biom file does not work, convert to JSON 
file = "Merged_otu_table.Nasal.Baseline.json.biom"
map = "180621_Map.OSA.Spain.Metadata.a4.BGW_edited.R.txt"

# Load the abundace table and mapping table 
abundance.table = import_biom(file, taxaPrefix=F)
mapping.table=sample_data(read.table(map, header=T, sep="\t", row.names=1))

lung.physeq=phyloseq(otu_table(abundance.table),tax_table(abundance.table), mapping.table)

#Give a colnames to separate different taxonomic levels
colnames(tax_table(lung.physeq))=c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "OTU")

# Load the tree file (use the unannotated.tree)
treefile = "97_otus_unannotated.tree"
tree.obj = import_qiime(treefilename = treefile) 

# Now merge the three separate phyloseq objects into a single object
otu.table = merge_phyloseq(lung.physeq, mapping.table, tree.obj)

rownames(sample_data(otu.table))
colnames(sample_data(otu.table))

# Remove taxa with 0 abundance
otu.table = subset_taxa(otu.table, rowSums(otu_table(otu.table)) != 0)

##If you want to nomalize OTU table before
## To normalize data you need to set a function
normalizeSample = function(x) {
    x/sum(x)
}
otu.relative.table = transformSampleCounts(otu.table, normalizeSample)

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))
otu.relative.table

rownames(sample_data(otu.table))
colnames(sample_data(otu.table))

Phylum.rel.table = tax_glom(otu.relative.table, taxrank = "Phylum")
Class.rel.table = tax_glom(otu.relative.table, taxrank = "Class")
Order.rel.table = tax_glom(otu.relative.table, taxrank = "Order")
Family.rel.table = tax_glom(otu.relative.table, taxrank = "Family")
Genus.rel.table = tax_glom(otu.relative.table, taxrank = "Genus")
OTU.rel.table = tax_glom(otu.relative.table, taxrank = "OTU")

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Zaragoza.RData")
load(file="Zaragoza.RData")

#################################################################
#################################################################
#################################################################

### To calculate beta diversity

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))
sample_data(otu.relative.table)$OSA_Cat

OSA.Zaragoza.relative.object = subset_samples(otu.relative.table, MSQ %in% c('56'))

### Calculate Beta diversity
### Calculate distance 
Unif.Zaragoza.phyloseq.dist = distance(OSA.Zaragoza.relative.object, method = "unifrac")

### estimate number of axes
UniF.Zaragoza.phyloseq.pco = dudi.pco(cailliez(Unif.Zaragoza.phyloseq.dist))

pdf(file="Zaragoza.allsamples.betadiversity.pdf", width=6, height=6)
s.class(UniF.Zaragoza.phyloseq.pco $li, interaction(sample_data(OSA.Zaragoza.relative.object)$OSA_Cat), col=c("forestgreen", "yellow3", "orange", "red"))
dev.off()

adonis(Unif.Zaragoza.phyloseq.dist~ OSA_Cat, data=data.frame(sample_data(OSA.Zaragoza.relative.object)), perm=1000)

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Zaragoza.RData")
load(file="Zaragoza.RData")

#################################################################
#################################################################
#################################################################

### All Samples

### Load Phyloseq
### Load all necessary packages 
library("phyloseq")
library("ggplot2")
library("ade4")
library("RColorBrewer")
library("ape")
library(plyr); library(dplyr)
library("plyr")
library("gplots")
library("RColorBrewer")
library("d3heatmap")
library("vegan")
library("igraph")
library("Heatplus")
library("graph")
library('ape')

### Install packages for heatmaps 
install.packages("gplots")
install.packages("RColorBrewer")
install.packages("d3heatmap")
install.packages("vegan")
install.packages("Heatplus")
install.packages("igraph")
source("https://bioconductor.org/biocLite.R")
biocLite("Heatplus")
source("https://bioconductor.org/biocLite.R")
biocLite("graph")
install.packages('gplots')
install.packages("ape")
install.packages("Heatplus")
install.packages("graphs")
install.packages("gplots")
install.packages("d3heatmap")

rownames(sample_data(Genus.rel.table))
colnames(sample_data(Genus.rel.table))

sample_data(Genus.rel.table)$SampleType
sample_data(Genus.rel.table)$Description


### Set of pruning we used:
### Here we can select for genera present in >2.5% relative abundance in 1.0% of the samples 
pdf(file="Allsamples.relative.abundance.pdf", width=20, height=8)
Genus.rel.OSA.table = genefilter_sample(OSA.Zaragoza.relative.object, filterfun_sample(function(x) x > 0.025), A = 0.01 * nsamples(OSA.Zaragoza.relative.object))
Genus.Rel.OSA.table.sequence = prune_taxa(Genus.rel.OSA.table, OSA.Zaragoza.relative.object)
plot_bar(Genus.Rel.OSA.table.sequence, fill="Genus")
dev.off()

rownames(sample_data(Genus.Rel.OSA.table.sequence))
colnames(sample_data(Genus.Rel.OSA.table.sequence))

#all data
data <- otu_table(OTU.rel.table)
#pruned to selected Genuses based on abundance
Genus.Data <-otu_table(Genus.Rel.OSA.table.sequence) 

#change colors to 0 = white, 1 = blue 
mypalette <- colorRampPalette(c("white","lightskyblue","navyblue"))(n=10)

### x10 into labRow for Genus.Rel.table.strep
### catglab for "genus"
### Here we are able to change the names for genuses that are labelled as "g__" --> Come back to this
x10 = prune_taxa(tail(names(sort(taxa_sums(Genus.Rel.OSA.table.sequence))), ntaxa(Genus.Rel.OSA.table.sequence)), Genus.Rel.OSA.table.sequence)
tax_table(x10)
# Add a new rank, Strain, with the Genus ids
tax_table(x10) <- cbind(tax_table(x10), Strain=taxa_names(x10))
# Define the ranks you want to include
myranks = c("Order", "Family", "Genus")
mylabels = apply(tax_table(x10)[, myranks], 1, paste, sep="", collapse="_")

# Add concatenated labels as a new rank after strain
tax_table(x10) <- cbind(tax_table(x10), catglab=mylabels)
# Check this out on a tree
plot_tree(x10, label.tips="catglab", nodelabf=nodeplotdefault, justify='jagged',ladderize="left", plot.margin=2.15, size="abundance")
plot_tree(x10, label.tips="catglab",ladderize="left", plot.margin=2.15, size="abundance")

###
###
###

breaks = seq(0,max(Genus.Data),length.out=1000)
gradient1 = colorpanel( sum( breaks[-0.01]<=0.1 ), "white", "lightskyblue")
gradient2 = colorpanel( sum( breaks[-1]>0.1 ), "lightskyblue", "navyblue" )
hm.colors = c(gradient1,gradient2)

distance1 = dist(Genus.Data)
cluster = hclust(distance1)

pdf(OSA.Eucleudian.Heatmap.AllSamples.pdf, height = 10, width = 20)
heatmap.2 (Genus.Data,
	density.info = "none",
	trace = "none",
	dendrogram = "both",
    labCol = sample_data(OSA.phyloseq.Genus.relative.object)$SampleIDUnique,
    labRow=tax_table(x10)[,"catglab"],
	cexCol = .35,
	cexRow = 1.0,
	symm=F,symkey=F,symbreaks=T, scale="none",
 	col=hm.colors,
 	breaks=breaks,
 	lhei = c(5,5),
 	margins=c(10,10),
)
dev.off()

#cluster Genuses(row)
GenusData.Bray.dist <-vegdist(Genus.Data, method = "bray")
Genus.Bray.clus <-hclust(GenusData.Bray.dist, "aver")

#cluster samples(Col)
Samples.Bray.dist = distance(Genus.Data, method="bray")
Samples.cluster.Bray = hclust(Samples.Bray.dist, "aver")
pdf("OSA.bray.Heatmap.AllSamples.pdf", height = 10, width = 20)
heatmap.2(Genus.Data, 
	density.info = "none",
	trace = "none",
	dendrogram = "both",
	Rowv = as.dendrogram(Genus.Bray.clus),
	Colv = as.dendrogram(Samples.cluster.Bray),
	labRow=tax_table(x10)[,"catglab"],
	cexCol = .35,
	cexRow = 1.0,
	labCol = sample_data(OSA.phyloseq.Genus.relative.object)$SampleIDUnique,
	col = hm.colors,
	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks=breaks,
 	lhei = c(5,5),
 	margins=c(10,10),
)
dev.off()

#Cluster taxa(rows) by phylogenetic origin. Because tree was not ultrametric, used chronos()
is.ultrametric(phy_tree(Genus.Rel.OSA.table.sequence))
Taxa_dend <- chronos(phy_tree(Genus.Rel.OSA.table.sequence))

#visualize tress to confirm they are the same 
plot(phy_tree(Genus.Rel.OSA.table.sequence), cex = .5) 
plot(Taxa_dend, cex = .5)

#Set the tree as a clustering object for import into heatmap
Taxa.clust <- as.hclust.phylo(Taxa_dend)

#Cluster Samples(col) by Unifrac
Samples.Unifrac.dist = UniFrac(Genus.Rel.OSA.table.sequence, weighted=TRUE)
Samples.cluster.Unifrac = hclust(Samples.Unifrac.dist, "aver")

pdf("OSA.unifrac.Heatmap.AllSamples.pdf", height = 10, width = 20)
heatmap.2(Genus.Data, 
	density.info = "none",
	trace = "none",
	Rowv = as.dendrogram(Taxa.clust),
	Colv = as.dendrogram(Samples.cluster.Unifrac),
	labRow=tax_table(x10)[,"catglab"],
	cexRow = .6,
	labCol = sample_data(Genus.Rel.OSA.table.sequence)$SampleIDUnique,
	cexCol = .35,
 	col = hm.colors,
	symm=F,symkey=F,symbreaks=T, scale="none",
 	breaks=breaks,
 	lhei = c(5,5),
 	margins=c(10,10),
)
dev.off()

###
###
###

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Zaragoza.RData")
load(file="Zaragoza.RData")

#################################################################
#################################################################
#################################################################

### Alpha diversity 

sessionInfo()
detach("package:igraph") 
sessionInfo()

### ALPHA DIVERSITY SHANNON 
### otu.relative.table

### All Samples Alpha diversity 
Shannon_diversity_all = diversity(otu_table(otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))

### Plot the diversity 
pdf(file="otu_diversity_Sample_type_allsamples.pdf", width=10, height=10)
boxplot(Shannon_diversity_all ~ sample_data(otu.relative.table)$OSA_Cat)
dev.off()
write.table(Shannon_diversity_all, file="Shannon_all.txt", sep="\t")

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Zaragoza.RData")
load(file="Zaragoza.RData")

#################################################################
#################################################################
#################################################################

### Calculate beta diversity
### Calculate distance 

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))
otu.relative.table

# Lymphlavagequartile_code (1, 0)
# IL8lavagequartile_code (1, 0)
# IL6lavagequartile_code (1, 0)

Lymphlavagequartile.Zaragoza.relative.object = subset_samples(OSA.Zaragoza.relative.object, Lymphlavagequartile_code %in% c('1', '0'))
IL8lavagequartile.Zaragoza.relative.object = subset_samples(OSA.Zaragoza.relative.object, IL8lavagequartile_code %in% c('1', '0'))
IL6lavagequartile.Zaragoza.relative.object = subset_samples(OSA.Zaragoza.relative.object, IL6lavagequartile_code %in% c('1', '0'))

Unif.Lymphlavagequartile_quartile.phyloseq.dist = distance(Lymphlavagequartile.Zaragoza.relative.object, method = "unifrac")
Unif.IL8lavagequartile.phyloseq.dist = distance(IL8lavagequartile.Zaragoza.relative.object, method = "unifrac")
Unif.IL6lavagequartile.phyloseq.dist = distance(IL6lavagequartile.Zaragoza.relative.object, method = "unifrac")

Unif.Lymphlavagequartile_quartile.phyloseq.pco = dudi.pco(cailliez(Unif.Lymphlavagequartile_quartile.phyloseq.dist), scannf = FALSE, nf = 3)
Unif.IL8lavagequartile.phyloseq.pco = dudi.pco(cailliez(Unif.IL8lavagequartile.phyloseq.dist), scannf = FALSE, nf = 3)
Unif.IL6lavagequartile.phyloseq.pco = dudi.pco(cailliez(Unif.IL6lavagequartile.phyloseq.dist), scannf = FALSE, nf = 3)

adonis(Unif.Lymphlavagequartile_quartile.phyloseq.dist ~ Lymphlavagequartile_code, data=data.frame(sample_data(Lymphlavagequartile.Zaragoza.relative.object)), perm=1000)
adonis(Unif.IL8lavagequartile.phyloseq.dist ~ IL8lavagequartile_code, data=data.frame(sample_data(IL8lavagequartile.Zaragoza.relative.object)), perm=1000)
adonis(Unif.IL6lavagequartile.phyloseq.dist ~ IL6lavagequartile_code, data=data.frame(sample_data(IL6lavagequartile.Zaragoza.relative.object)), perm=1000)

all.ordinate.Zaragoza.Lymphlavagequartile.UniF <- ordinate(Lymphlavagequartile.Zaragoza.relative.object, method="PCoA", distance="unifrac")
all.ordinate.Zaragoza.IL8lavagequartile.UniF <- ordinate(IL8lavagequartile.Zaragoza.relative.object, method="PCoA", distance="unifrac")
all.ordinate.Zaragoza.IL6lavagequartile.UniF <- ordinate(IL6lavagequartile.Zaragoza.relative.object, method="PCoA", distance="unifrac")

pdf(file="UniFrac.ordination.Zaragoza.Lymphlavagequartile.OSA.pdf", width=6, height=6) 
plot_ordination(Lymphlavagequartile.Zaragoza.relative.object, all.ordinate.Zaragoza.Lymphlavagequartile.UniF, color = "OSA_Cat", shape = "OSA_Cat")
dev.off()
pdf(file="UniFrac.ordination.Zaragoza.IL8lavagequartile.OSA.pdf", width=6, height=6) 
plot_ordination(IL8lavagequartile.Zaragoza.relative.object, all.ordinate.Zaragoza.IL8lavagequartile.UniF, color = "OSA_Cat", shape = "OSA_Cat")
dev.off()
pdf(file="UniFrac.ordination.Zaragoza.IL6lavagequartile.OSA.pdf", width=6, height=6) 
plot_ordination(IL6lavagequartile.Zaragoza.relative.object, all.ordinate.Zaragoza.IL6lavagequartile.UniF, color = "OSA_Cat", shape = "OSA_Cat")
dev.off()

pdf(file="UniFrac.Zaragoza.Lymphlavagequartile.OSA.pdf", width=6, height=6) 
s.class(Unif.Lymphlavagequartile_quartile.phyloseq.pco $li, interaction(sample_data(Lymphlavagequartile.Zaragoza.relative.object)$Lymphlavagequartile), col=c("red", "forestgreen"), label=c('Lymph high', 'Lymph low'))
dev.off()
pdf(file="UniFrac.Zaragoza.IL8lavagequartile.OSA.pdf", width=6, height=6) 
s.class(Unif.IL8lavagequartile.phyloseq.pco $li, interaction(sample_data(IL8lavagequartile.Zaragoza.relative.object)$IL8lavagequartile), col=c("red", "forestgreen"), label=c('IL8 high', 'IL8 low'))
dev.off()
pdf(file="UniFrac.Zaragoza.IL6lavagequartile.OSA.pdf", width=6, height=6) 
s.class(Unif.IL6lavagequartile.phyloseq.pco $li, interaction(sample_data(IL6lavagequartile.Zaragoza.relative.object)$IL6lavagequartile), col=c("red", "forestgreen"), label=c('IL6 high', 'IL6 low'))
dev.off()

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Zaragoza.RData")
load(file="Zaragoza.RData")

#################################################################
#################################################################
#################################################################