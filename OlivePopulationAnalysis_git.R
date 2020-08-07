#########################################################################################
## SET UP
#########################################################################################
## Set up the working directory
setwd("")

## Call the libraries
library(vcfR)
library(poppr)
library(ape)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)

library(fields)
library(mapplots)
library(LEA)

library(SNPRelate)
library(dartR)

## Load some LEA functions
source("http://membres-timc.imag.fr/Olivier.Francois/Conversion.R")
source("http://membres-timc.imag.fr/Olivier.Francois/POPSutilities.R")

## Keep the default margins before we change them
default_par_mar = par("mar")
default_par_oma = par("oma")

#########################################################################################
## LOAD THE SAMPLES FOR THE POPULATION STRUCTURE ANALYSIS
#########################################################################################

## 1.1- Load the sample

OleaPanel = read.vcfR("C:/file/location/filename.vcf")

## Generate different data formats and upload the population information
GL_OleaPanel = vcfR2genlight(OleaPanel)

pop.data = read.csv("C:/file/location/filename.csv", sep = ",", header = TRUE)
colnames(pop.data) = c("AccessID", "Origin", "Group")
all(colnames(OleaPanel@gt)[-1] == pop.data$AccessID)
pop(GL_OleaPanel) = pop.data$Origin



gl2faststructure(GL_OleaPanel, outfile = "OleurDiversityPanel.Variants20200412.fs.str")

#########################################################################################
## STRUCTURE ANALYSIS WITH FASTSTRUCTURE
#########################################################################################

## 2.1- Create the geno File
struct2geno(file = "OleurDiversityPanel.Variants20200412.fs.str", TESS = FALSE, diploid = TRUE, FORMAT = 2, extra.row = 0, extra.col = 0, output = "./OleurDiversityPanel.Variants20200412.geno")

## 2.2- Run the Structure analysis with number of clusters (K) from 1 to 20
obj_str = snmf("./OleurDiversityPanel.Variants20200412.geno", K = 1:20, ploidy = 2, entropy = T, CPU = 2, project = "new")

## 2.3- Check for the best K
plot(obj_str, col = "blue4", cex = 1.4, pch = 19)

## Results: It looks like two possible best K are 2 and 5 (save the plot for supplementals)

## 2.4- Create two barplots, at K=2 and K=5
qmatrix_k3 = Q(obj_str, K =2)
par(oma=c(5,0,1,0), mar=c(5,5,0,0))
barplot(t(qmatrix_k3), col = c("lightblue", "violet", "blue"), border = NA, space = 0, ylab = "Admixture coefficients", names.arg = GL_OleaPanel$ind.names, las=2, cex.names = 0.7)

qmatrix_k5 = Q(obj_str, K =5)
par(oma=c(5,0,1,0), mar=c(5,5,0,0))
barplot(t(qmatrix_k5), col = c("lightblue", "violet", "blue", "red", "pink"), border = NA, space = 0, ylab = "Admixture coefficients", names.arg = GL_OleaPanel$ind.names, las=2, cex.names = 0.7)

#########################################################################################
## DISTANCE TREE
#########################################################################################
## 3.1- Generate a distance tree

Olea_tree = aboot(GL_OleaPanel, tree = "nj", distance = bitwise.dist, sample = 100, showtree = FALSE, cutoff = 50, quiet = T)

rainbow(n = nPop(GL_OleaPanel))
cols = rainbow(n = nPop(GL_OleaPanel))

par(oma=c(3,2,1,1), mar=c(0,0,0,0))
plot.phylo(Olea_tree, cex = 0.8, font = 2, adj = 0, tip.color =  cols[pop(GL_OleaPanel)])
nodelabels(Olea_tree$node.label, adj = c(1.3, -0.5), frame = "n", cex = 0.8,font = 3, xpd = TRUE)
legend('bottomright', legend = levels(pop(GL_OleaPanel)), fill = cols, border = FALSE, bty = "n", cex = 0.7)
axis(side = 1)
title(xlab = "Genetic distance (proportion of loci that are different)")


#########################################################################################
## PCA ANALYSIS
#########################################################################################
## 4.1- Run PCA analysis

Olea_PCA = glPca(GL_OleaPanel_Rm, nf=3)

## 4.2-Plot the PCA

## Set up new margins
par(oma=c(3,3,3,3), mar=c(1,1,1,1))

## Plot the barplot for the eigenvalues
barplot(100*Olea_PCA$eig/sum(Olea_PCA$eig), col = heat.colors(50), main="PCA Eigenvalues")
title(ylab="Percent of variance\nexplained", line = 2)
title(xlab="Eigenvalues", line = 1)
Olea_PCA_scores = as.data.frame(Olea_PCA$scores)
Olea_PCA_scores$pop = pop(GL_OleaPanel_Rm)

## Plot the PCA
p <- ggplot(Olea_PCA_scores, aes(x=PC1, y=PC2, colour=pop))
p <- p + geom_point(size=2)

## To show the labels we will use the geom_text_repel function
p <- p + geom_text_repel(aes(label=gsub("OE._..._", "", rownames(Olea_PCA_scores))), show.legend = FALSE)
p <- p + scale_color_manual(values = cols)
p <- p + geom_hline(yintercept = 0)
p <- p + geom_vline(xintercept = 0)
p <- p + theme_bw()
p

#########################################################################################
## DAPC ANALYSIS
#########################################################################################

## 5.0- Call the libraries
library(reshape2)

## 5.1- Run the DAPC analysis with 3 components and groups based 
Olea_DAPC <- dapc(GL_OleaPanel_Rm, n.pca = 3, n.da = 2)

## 5.2- Plot components
scatter(Olea_DAPC, col = cols, cex = 2, legend = TRUE, clabel = F, posi.leg = "bottomleft", scree.pca = TRUE, posi.pca = "topleft", cleg = 0.75)
compoplot(Olea_DAPC,col = cols, posi = 'top')

## 5.3- Generate the barplot with populations assigned to each sector
dapc.results <- as.data.frame(Olea_DAPC$posterior)
dapc.results$pop <- pop(GL_OleaPanel_Rm)
dapc.results$indNames <- rownames(dapc.results)

dapc.results <- melt(dapc.results)
colnames(dapc.results) <- c("Original_Pop","Sample","Assigned_Pop","Posterior_membership_probability")

p <- ggplot(dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Pop))
p <- p + geom_bar(stat='identity')
p <- p + scale_fill_manual(values = cols)
p <- p + facet_grid(~Original_Pop, scales = "free")
p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p

## The DAPC can be used also to infer the number of groups and then plot them
## in a similar fashion than the Structure. To do it:

## 5.4- Check the number of groups and their variation to find the optimal K

## Use as maximum K=10 and create a matrix with the Kstat results of find clusters 
maxK = 10
myMat = matrix(nrow=10, ncol=maxK)
colnames(myMat) = 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp = find.clusters(GL_OleaPanel_Rm, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = maxK)
  myMat[i,] <- grp$Kstat
}

## Generate a data.frame for the plot
my_df = melt(myMat)
colnames(my_df)[1:3] = c("Group", "K", "BIC")
my_df$K = as.factor(my_df$K)
head(my_df)

## Plot the BIC values associated with find.clusters
p1 = ggplot(my_df, aes(x = K, y = BIC))
p1 = p1 + geom_boxplot()
p1 = p1 + theme_bw()
p1 = p1 + xlab("Number of groups (K)")
p1

## Results: It looks like the optimal number of groups is 2 (like in the Structure analysis)

## 5.5- Create the accession group assignments based in the LD 

## We will assay from 2 to 5 groups
my_k = 2:5

grp_l = vector(mode = "list", length = length(my_k))
dapc_l = vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  set.seed(9)
  grp_l[[i]] = find.clusters(GL_OleaPanel_Rm, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] = dapc(GL_OleaPanel_Rm, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
}

my_df = as.data.frame(dapc_l[[ length(dapc_l) ]]$ind.coord)
my_df$Group = dapc_l[[ length(dapc_l) ]]$grp
head(my_df)

my_pal = RColorBrewer::brewer.pal(n=8, name = "Dark2")

p2 = ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group))
p2 = p2 + geom_point(size = 4, shape = 21)
p2 = p2 + theme_bw()
p2 = p2 + scale_color_manual(values=c(my_pal))
p2 = p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2

## 5.6- Create the accession group assignments based in the posterior probability
tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- pop(GL_OleaPanel_Rm)
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- pop(GL_OleaPanel_Rm)
  
  my_df <- rbind(my_df, tmp)
}

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(my_pal))
p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
p3

#########################################################################################
## BASIC STATS FOR GROUPS (based in the groups produced by DAPC)
#########################################################################################

## Once the different individuals have assigned to different populations, last step will be to redo the population 
## assignment and calculate some population genetic parameters.

my_dfK2 = my_df[my_df$K == 3,]
K3groups = my_dfK2[my_dfK2$Posterior == 1,]
GL_OleaPanel_RmCp = GL_OleaPanel_Rm
pop(GL_OleaPanel_RmCp) = K3groups$Group

GL_OleaPanel_RmCp_BStats = gl.basic.stats(GL_OleaPanel_RmCp)


#########################################################################################
## BASIC STATS USING POPGENOME
#########################################################################################

### It is necessary to break the VCF into sequences (for readData). To do it it is necessary to use bcftools
### in Linux:
### 1- Compress the file: "bgzip -c myfile.vcf > myfile.vcf.gz"
### 2- Index the compressed file "tabix -p vcf myfile.vcf.gz"
### 3- Create a script to divide the VCF file per sequences
###    grep -v "#" myfile.vcf | cut -f1 | sort -u | awk 'BEGIN{ print "#!/bin/bash\n\nmkdir vcf_by_seq\n"}{ print "bcftools view myfile.vcf.gz "$1" > vcf_by_seq/"$1}' > RunGetVcfBySeq.sh
### 4- Executate the script   

library(PopGenome)
## 6.1- Load the data
PG_OleaPanel = readData("C:/file/location/filename", format = "VCF")


## Check the total number of SNPs
PG_OleaPanel_sumdata = as.data.frame(get.sum.data(PG_OleaPanel))
sum(PG_OleaPanel_sumdata$n.biallelic.sites)

## Add the populations
K3group1list = K3groups[K3groups$Group == 1,1]
K3group2list = K3groups[K3groups$Group == 2,1]
K3group3list = K3groups[K3groups$Group == 3,1]
PG_OleaPanel = set.populations(PG_OleaPanel, list(K3group1list, K3group2list, K3group3list), diploid = TRUE)

## 6.2- Estimate the different parameters
PG_OleaPanel_concatenated = concatenate.regions(PG_OleaPanel)
PG_OleaPanel_concatenated = neutrality.stats(PG_OleaPanel_concatenated, FAST=TRUE)
get.neutrality(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = diversity.stats(PG_OleaPanel_concatenated)
get.diversity(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = F_ST.stats(PG_OleaPanel_concatenated)
get.F_ST(PG_OleaPanel_concatenated)
PG_OleaPanel_concatenated = detail.stats(PG_OleaPanel_concatenated, site.spectrum=TRUE, site.FST=TRUE) 
PG_OleaPanel_results = get.detail(PG_OleaPanel_concatenated) 

## Get the segregating sites
PG_OleaPanel_concatenated@n.segregating.sites

## Get nucleotide diversity, as Pi, for each population
PG_OleaPanel_concatenated@Pi

## Get the Tajima D for each of the populations
PG_OleaPanel_concatenated@Tajima.D

## Get the Waterson Theta for each of the populations
PG_OleaPanel_concatenated@theta_Watterson

> 
