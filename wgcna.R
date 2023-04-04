library(WGCNA)
if (!require("BiocManager" , quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
if (!require("BiocManager" , quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GO.db")
if (!require("BiocManager" , quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("lib")

install.packages("BiocManager")
BiocManager::install("WGCNA")

library(WGCNA)
install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"));

library(WGCNA)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GO.db")

if (!require("BiocManager" , quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
library(DESeq2)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("tidyverse")
library(GEOquery)
source("http://www.bioconductor.org/biocLite.R")
biocLite(c("GEOquery"))
## try http:// if https:// URLs are not supported
source("https://bioconductor.org/biocLite.R")
biocLite("GEOquery")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CorLevelPlot")


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("gridExtra")

library(WGCNA)
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("GO.db")

library(WGCNA)


library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)

#install.packages("CorLevelPlot")
#devtools::install_github("kevinblighe/CorLevelPlot")
??WGCNAnThreads
allowWGCNAThreads()
data <- read.delim('/Users/sabinamahnesaei/Downloads/data.txt', header = T)
geo_id <- "GSE152418"
gse <- getGEO(geo_id, GSEMatrix = TRUE)
phenoData <- pData(phenoData(gse[[1]]))

#SHOW= head
head(phenoData)


#bazaro kochektar kardim
phenoData <- phenoData [,c(1,2, 46:50)]




data [1:10, 1:10]


#change the sample into gene accession number
data<- data %>%
  gather(key = 'samples', value = 'counts', -ENSEMBLID) %>%
  mutate(samples = gsub ('\\.', '-' , samples)) %>%
  inner_join(., phenoData, by = c('samples'= 'title')) %>%
  #baraye inke ye seri az bakhsh haro joda koni az code zir estefadeh kon
  
  select(1,3,4) %>%
  
  #change the shape
  spread(key = 'geo_accession', value = 'counts') %>%
  column_to_rownames(var = 'ENSEMBLID')

gsg <- goodSamplesGenes(t (data))
  summary(gsg)
  gsg$allOK
table (gsg$goodGenes)  
table(gsg$goodSamples)
data <- data [gsg$goodGenes == TRUE ,]
htree <- hclust(dist(t(data)), method = "average")
plot(htree)
#another method
pca <-prcomp(t(data))
pca.dat <- pca$sdev
pca.var <-pca$sdev ^ 2
pca.var.percent <- round(pca.var/sum(pca.var)*100 , digits = 2)


pca.dat <- as.data.frame(pca.dat)

ggplot(pca.dat, aes(PC1, PC2)) +
  geom_point() +
  geom_text(label = rownames(pca.dat)) +
  labs(x = paste0('pc1: ', pca.var.percent[1], '%'),
       y = paste0('pc2: ', pca.var.percent[2], '%'))
 
library(plotly)


#exclude outlier samples

samples.to.be.excluded <- c('GSM1615000', 'GSM4614993' , 'GSM4614995')
data.subset <- data[,!(colnames(data) %in% samples.to.be.excluded)]


# 3.normalization
#phenoData baraye expression estefadeh mishe
#create a deseq2 
#DESeq2 = The DESeq2 package in R is used for differential gene expression analysis
#based on the negative binomial distribution and use high throuput from sequences like RNA seq

colData <-phenoData %>%
  filter(!row.names(.) %in% samples.to.be.excluded)
#fixing colnames in coldata
names(colData)


#gsub = baraye remove space va charecter estefadeh mishe
names(colData) <- gsub (':ch1', '', names(colData) )
names(colData) <-gsub ('\\s', '_', names(colData))
#making the rownames and column names in indentical

all(rownames(colData) %in% colnames(data.subset))
all(rownames(colData) == colnames(data.subset))


# create dds
dds <- DESeqDataSetFromMatrix( countData = data.subset,
                               colData = colData,
                               design = ~ 1) #not specify models
##remove all genes with counts < 15 in more than 75% of samples
## suggested by wgcna in rna seq FAQ

dds75 <- dds[rowSums(counts(dds) >= 15) >=24,]
nrow(dds75) #13360

#perform various stabilazition
dds_norm <- vst(dds75)
#get normalization counts
norm.counts <- assay(dds_norm) %>%
  t()
  #head()


# 4.Network construction
# choose a set of soft thershoulding powers



power <- c(c(1:10),  seq(from =20, to =50, by=2))
power
#call the network topology analysis function
ls()

sft <-pickSoftThreshold( norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)
# hatman ghesmat networktype " " ino bezar
sft.data <- sft$fitIndices


#visualisation to pick power (max to min)

a1 <- ggplot(sft.data, aes(power,SFT.R.sq, label = power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color= 'red') +
  labs(x= 'Power', y = 'sclale free topology model fit, signed R^2') +
  theme_classic()



a2 <- ggplot(sft.data , aes(power,mean.k., label= power )) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs (x= 'power', y ='Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow(2)) 

# convert matrix to numeric


norm.counts [] <- sapply(norm.counts, as.numeric)
soft_power <- 20  
temp_cor <- cor
cor<- WGCNA::cor


#memory estimate w.r.t blocksize
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 5000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)
#hatmn to ghesmat beshin ye seri maghalat bekhon bebin kodom javab khobi azash migiri


cor <- temp_cor

# 5. module eigengens

module_Eigengenes <- bwnet$MEs
#print out the perview
head(module_Eigengenes)



#get number of gene for each modules


table(bwnet$colors)


#plot the dendogram and the module colors befor and after merging underneeth


plotDendroAndColors(bwnet$dendrograms[[1]], cbind(bwnet$unmergedColors, bwnet$colors),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    
                    addGuide = TRUE,
                    hang = 0.03,
                    guideHang = 0.05)
