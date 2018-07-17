### Scripts used in Microbiome analysis

update.packages()
getwd()
setwd('### working directory should go here')

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
biom convert -i /Volumes/homes/ben/Projects/SNORE-3/Merge/Closed/R/otu_table_new.filter.merge.biom -o /Volumes/homes/ben/Projects/SNORE-3/Merge/Closed/R/otu_table_new.filter.merge.txt --header-key "taxonomy" --to-tsv

biom convert -i /Volumes/homes/ben/Projects/SNORE-3/Merge/Closed/R/otu_table_new.filter.merge.txt -o /Volumes/homes/ben/Projects/SNORE-3/Merge/Closed/R/otu_table_new.filter.merge.json.biom --to-json --table-type="OTU table" --process-obs-metadata taxonomy
### QIIME convert to JSON ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### Load the files needed, if biom file does not work, convert to JSON 
file = "otu_table_new.filter.merge.json.biom"
map = "Merged.Snore.Map.ba.Final.BGW.upload.txt"

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
save.image(file="Snore.code.RData")
load(file="Snore.code.RData")

#################################################################
#################################################################
#################################################################

### To calculate beta diversity
### Prints out the available distance methods under UniFrac

rownames(sample_data(otu.relative.table))
colnames(sample_data(otu.relative.table))

OSA.phyloseq.relative.object = subset_samples(otu.relative.table, OSA_DX_4 %in% c('No', 'Mild', 'Moderate', 'Severe'))

### Calculate beta diversity
### Calculate distance 
Unif.OSA.phyloseq.dist = distance(OSA.phyloseq.relative.object, method = "unifrac")

### estimate number of axes
UniF.OSA.phyloseq.pco = dudi.pco(cailliez(Unif.OSA.phyloseq.dist))

pdf(file="OSA.category.expsosure.pdf", width=6, height=6)
s.class(UniF.OSA.phyloseq.pco $li, interaction(sample_data(OSA.phyloseq.relative.object)$OSA_DX_4), col=c("forestgreen", "yellow3", "orange", "red"))
dev.off()

adonis(Unif.OSA.phyloseq.dist ~ OSA_DX_4, data=data.frame(sample_data(OSA.phyloseq.relative.object)), perm=1000)

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Snore.code.RData")
load(file="Snore.code.RData")

#################################################################
#################################################################
#################################################################

### All Samples
### Heatmaps

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

OSA.phyloseq.Genus.relative.object = subset_samples(Genus.rel.table, OSA_DX_4 %in% c('No', 'Mild', 'Moderate', 'Severe'))

### Set of pruning we used:
### Here we can select for genera present in >2.5% relative abundance in 1.0% of the samples
pdf(file="Allsamples.relative.abundance.pdf", width=20, height=8)
Genus.rel.OSA.table = genefilter_sample(OSA.phyloseq.Genus.relative.object, filterfun_sample(function(x) x > 0.025), A = 0.01 * nsamples(OSA.phyloseq.Genus.relative.object))
Genus.Rel.OSA.table.sequence = prune_taxa(Genus.rel.OSA.table, Genus.rel.table)
plot_bar(Genus.Rel.OSA.table.sequence, fill="Genus")
dev.off()

rownames(sample_data(Genus.Rel.OSA.table.sequence))
colnames(sample_data(Genus.Rel.OSA.table.sequence))

sample_data(Genus.Rel.OSA.table.sequence)$Site
sample_data(Genus.Rel.OSA.table.sequence)$Site_Code

# All data
data <- otu_table(OTU.rel.table)
#pruned to selected Genuses based on abundance
Genus.Data <-otu_table(Genus.Rel.OSA.table.sequence) 

# Change colors to 0 = white, 1 = blue 
mypalette <- colorRampPalette(c("white","lightskyblue","navyblue"))(n=10)

# Add Col Side Colors by Sample Location
# create vector to lable by Topograph Code
SampleTypeVector = sample_data(Genus.Rel.OSA.table.sequence)$Site_Code

# Duplicate to create a color vector and replace value w/ color 
# Colorvector.strep can only be numbers
Colorvector <- SampleTypeVector

### NYU
Colorvector <- replace(Colorvector, which (Colorvector == "1"), "purple")
### Rutgers 
Colorvector <- replace(Colorvector, which (Colorvector == "2"), "goldenrod")

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

pdf("OSA.Eucleudian.Heatmap.AllSamples.pdf", height = 10, width = 20)
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
    ColSideColors= Colorvector,
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
	ColSideColors = Colorvector,
	main = "Heatmap of Bray-Curtis Distance",
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
	ColSideColors = Colorvector,
)
dev.off()

###
###
###

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Snore.code.RData")
load(file="Snore.code.RData")

#################################################################
#################################################################
#################################################################

### ALPHA DIVERSITY SHANNON 

### All Samples Alpha diversity 
Shannon_diversity_all = diversity(otu_table(otu.relative.table), index = "shannon", MARGIN = 2, base = exp(1))

### Plot the diversity 
pdf(file="otu_diversity_Sample_type_allsamples.pdf", width=10, height=10)
boxplot(Shannon_diversity_all ~ sample_data(otu.relative.table)$OSA_DX_4)
dev.off()
write.table(Shannon_diversity_all, file="Shannon_all.txt", sep="\t")

OSA.phyloseq.relative.object = subset_samples(otu.relative.table, OSA_DX_4 %in% c('No', 'Mild', 'Moderate', 'Severe'))
Total.Neut.mL_quartile.relative.object = subset_samples(otu.relative.table, Total.Neut.mL_quartile_code %in% c('1', '0'))
IL_8_quartile.relative.object = subset_samples(otu.relative.table, IL_8_quartile_code %in% c('1', '0'))
IL_6_quartile.relative.object = subset_samples(otu.relative.table, IL_6_quartile_code %in% c('1', '0'))

Shannon_diversity_OSA = diversity(otu_table(OSA.phyloseq.relative.object), index = "shannon", MARGIN = 2, base = exp(1))
Shannon_diversity_Neut = diversity(otu_table(Total.Neut.mL_quartile.relative.object), index = "shannon", MARGIN = 2, base = exp(1))
Shannon_diversity_IL8 = diversity(otu_table(IL_8_quartile.relative.object), index = "shannon", MARGIN = 2, base = exp(1))
Shannon_diversity_IL6 = diversity(otu_table(IL_6_quartile.relative.object), index = "shannon", MARGIN = 2, base = exp(1))

newSTorder = c("No", "Mild", "Moderate", "Severe")
sample_data(OSA.phyloseq.relative.object)$OSA_DX_4 <- as.character(sample_data(OSA.phyloseq.relative.object)$OSA_DX_4)
sample_data(OSA.phyloseq.relative.object)$OSA_DX_4 <- factor(sample_data(OSA.phyloseq.relative.object)$OSA_DX_4, levels=newSTorder)

pdf(file="all_samples.alpha.OSA.pdf", width=10, height=6)
boxplot(Shannon_diversity_OSA ~ sample_data(OSA.phyloseq.relative.object)$OSA_DX_4, col=c("forestgreen", "gold3", "Orange3","red4"), outline=FALSE)
stripchart(Shannon_diversity_OSA ~ sample_data(OSA.phyloseq.relative.object)$OSA_DX_4, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()

pdf(file="Total.Neut.mL.alpha.OSA.pdf", width=6, height=6)
boxplot(Shannon_diversity_Neut ~ sample_data(Total.Neut.mL_quartile.relative.object)$Total.Neut.mL_quartile, col=c("red4", "forestgreen"), outline=FALSE)
stripchart(Shannon_diversity_Neut ~ sample_data(Total.Neut.mL_quartile.relative.object)$Total.Neut.mL_quartile, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()

pdf(file="IL_8_quartile.alpha.OSA.pdf", width=6, height=6)
boxplot(Shannon_diversity_IL8 ~ sample_data(IL_8_quartile.relative.object)$IL_8_quartile, col=c("red4", "forestgreen"), outline=FALSE)
stripchart(Shannon_diversity_IL8 ~ sample_data(IL_8_quartile.relative.object)$IL_8_quartile, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()

pdf(file="IL_6_quartile.alpha.OSA.pdf", width=6, height=6)
boxplot(Shannon_diversity_IL6 ~ sample_data(IL_6_quartile.relative.object)$IL_6_quartile, col=c("red4", "forestgreen"), outline=FALSE)
stripchart(Shannon_diversity_IL6 ~ sample_data(IL_6_quartile.relative.object)$IL_6_quartile, vertical=TRUE, add=TRUE, method="jitter", pch=20, col="black")
dev.off()

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Snore.code.RData")
load(file="Snore.code.RData")

#################################################################
#################################################################
#################################################################

Total.Neut.mL_quartile.relative.object = subset_samples(otu.relative.table, Total.Neut.mL_quartile_code %in% c('1', '0'))
IL_8_quartile.relative.object = subset_samples(otu.relative.table, IL_8_quartile_code %in% c('1', '0'))
IL_6_quartile.relative.object = subset_samples(otu.relative.table, IL_6_quartile_code %in% c('1', '0'))

Unif.Neut.mL_quartile.phyloseq.dist = distance(Total.Neut.mL_quartile.relative.object, method = "unifrac")
Unif.IL_8_quartile.phyloseq.dist = distance(IL_8_quartile.relative.object, method = "unifrac")
Unif.IL_6_quartile.phyloseq.dist = distance(IL_6_quartile.relative.object, method = "unifrac")

UniF.Neut.mL_quartile.phyloseq.pco = dudi.pco(cailliez(Unif.Neut.mL_quartile.phyloseq.dist), scannf = FALSE, nf = 3)
UniF.IL_8_quartile.phyloseq.pco = dudi.pco(cailliez(Unif.IL_8_quartile.phyloseq.dist), scannf = FALSE, nf = 3)
UniF.IL_6_quartile.phyloseq.pco = dudi.pco(cailliez(Unif.IL_6_quartile.phyloseq.dist), scannf = FALSE, nf = 3)

adonis(Unif.Neut.mL_quartile.phyloseq.dist ~ Total.Neut.mL_quartile_code, data=data.frame(sample_data(Total.Neut.mL_quartile.relative.object)), perm=1000)
adonis(Unif.IL_8_quartile.phyloseq.dist ~ IL_8_quartile_code, data=data.frame(sample_data(IL_8_quartile.relative.object)), perm=1000)
adonis(Unif.IL_6_quartile.phyloseq.dist ~ IL_6_quartile_code, data=data.frame(sample_data(IL_6_quartile.relative.object)), perm=1000)

pdf(file="Total.Neut.mL.pdf", width=6, height=6)
s.class(UniF.Neut.mL_quartile.phyloseq.pco$li, interaction(sample_data(Total.Neut.mL_quartile.relative.object)$Total.Neut.mL_quartile), col=c("red", "forestgreen"))
dev.off()
pdf(file="IL_8_quartile.pdf", width=6, height=6)
s.class(UniF.IL_8_quartile.phyloseq.pco$li, interaction(sample_data(IL_8_quartile.relative.object)$IL_8_quartile), col=c("red", "forestgreen"))
dev.off()
pdf(file="IL_6_quartile.pdf", width=6, height=6)
s.class(UniF.IL_6_quartile.phyloseq.pco$li, interaction(sample_data(IL_6_quartile.relative.object)$IL_6_quartile), col=c("red", "forestgreen"))
dev.off()

#GET AXIS
all.ordinate.Neut.UniF <- ordinate(Total.Neut.mL_quartile.relative.object, method="PCoA", distance="unifrac")
all.ordinate.IL8.UniF <- ordinate(IL_8_quartile.relative.object, method="PCoA", distance="unifrac")
all.ordinate.IL6.UniF <- ordinate(IL_6_quartile.relative.object, method="PCoA", distance="unifrac")

# Scatter plot
pdf(file="ordination.Neut.exposure.pdf") 
plot_ordination(Total.Neut.mL_quartile.relative.object, all.ordinate.Neut.UniF, color = "OSA_DX_4", shape = "OSA_DX_4")
dev.off()
pdf(file="ordination.IL8.exposure.pdf") 
plot_ordination(IL_8_quartile.relative.object, all.ordinate.IL8.UniF, color = "OSA_DX_4", shape = "OSA_DX_4")
dev.off()
pdf(file="ordination.IL6.exposure.pdf") 
plot_ordination(IL_6_quartile.relative.object, all.ordinate.IL6.UniF, color = "OSA_DX_4", shape = "OSA_DX_4")
dev.off()

#################################################################
#################################################################
#################################################################

otu.relative.table
save.image(file="Snore.code.RData")
load(file="Snore.code.RData")

#################################################################
#################################################################
#################################################################
