# activating all the packages needed for the code

#provides functions for handling hdf5 file formats (kallisto outputs bootstraps in this format)
library(rhdf5) 

# provides access to Hadley Wickham's collection of R packages for data science, which we will use throughout the course
library(tidyverse) 

# package for getting Kallisto results into R
library(tximport) 

#helps deal with ensembl
library(ensembldb)

# for the annotations
library(biomaRt)

# well known package for differential expression analysis, but we only use for the DGEList object and for normalization methods
library(edgeR) 

# let's us easily calculate stats on rows or columns of a data matrix
library(matrixStats)

# allows you to combine multiple plots in one figure
library(cowplot) 

library(ggplot2)

# for making interactive plots
library(plotly)

# A layered 'grammar of tables' - think ggplot, but for tables
library(gt)

library(reshape2)

library(limma)

#for heatmaps
library(gplots)

#interactive and searchable tables of our GSEA results
library(DT)

#functions and methods for Gene Set Enrichment Analysis
library(GSEABase)
        
#base functions for bioconductor; required by GSEABase
library(Biobase)

#Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment 
# across samples.
library(GSVA) 

#tools for accessing the GO enrichment results using g:Profiler web resources
library(gprofiler2) 

# provides a suite of tools for functional enrichment analysis
library(clusterProfiler) 

# access to msigdb collections directly within R
library(msigdbr) 

# great for making the standard GSEA enrichment plots
library(enrichplot)

# Quantitative Set Analysis for Gene Expression
library(qusage) 

library(heatmaply)
library(dplyr)
library(patchwork)

# library(biocLite) I think I need this but I'm not completely sure for what

################################################################################################################################

# reading in the study design
targets <- read_tsv("studydesign.txt")

# set file paths to mapped data
path <- file.path(targets$sra_accession, "abundance.tsv")

# checking if it worked
all(file.exists(path))

# getting the ensembl annotations 
rice_anno <- useMart(biomart = "plants_mart", dataset = "osativa_eg_gene", host ="https://plants.ensembl.org")
rice_attributes <- listAttributes(rice_anno)

# retrieving user-specific attributes, we want "ensembl_transcript_id" and "external_gene_name" and we're extracting
# it from our mart rice.anno where we saved our annotations and convert it into a tibble (dataframe)
rice <- getBM(attributes=c('ensembl_transcript_id',
                              'external_gene_name'),
                 mart = rice_anno)

rice <- as_tibble(rice)

#renaming the columns retrieved from bioMart
rice <- dplyr::rename(rice, target_id = ensembl_transcript_id, 
                         gene_name = external_gene_name)

#importing the Kallisto output into R
Txi_gene <- tximport(path, 
                     type = "kallisto", 
                     tx2gene = rice, 
                     txOut = FALSE, 
                     countsFromAbundance = "lengthScaledTPM",
                     ignoreTxVersion = TRUE)

# Transcripts per million (TPM) is a measurement of the proportion of transcripts in your pool of RNA.
myTPM <- Txi_gene$abundance

# “Counts” usually refers to the number of reads that align to a particular feature
myCounts <- Txi_gene$counts

#capture sample labels from the study design file
sampleLabels <- targets$sample

# Generating summary stats for the data
# adds standard deviation, average and median of every row to myTPM matrix
myTPM_stats <- transform(myTPM, 
                         SD=rowSds(myTPM), 
                         AVG=rowMeans(myTPM),
                         MED=rowMedians(myTPM))

# Making a DGElist from counts
# Differential gene expression: A gene is declared differentially expressed if an observed difference or change in read counts 
# or expression levels between two experimental conditions is statistically significant.
DGE_list <- DGEList(myCounts)

# getting counts per million with the 'cpm' function (Counts per million (CPM) mapped reads are counts scaled by the number of 
# fragments you sequenced (N) times one million)
cpm_of_DGE <- cpm(DGE_list) 

# this takes the log with base two of the given number in cpm
log2_cpm_of_DGE <- cpm(DGE_list, log=TRUE)

# 'coerce' your data matrix to a dataframe so that you can use tidyverse tools on it
log2_cpm_of_DGE_df <- as_tibble(log2_cpm_of_DGE, rownames = "geneID")

# add your sample names to this dataframe (we lost these when we read our data in with tximport)
colnames(log2_cpm_of_DGE_df) <- c("geneID", sampleLabels)

# use the tidy package to 'pivot' your dataframe (from wide to long)
log2_cpm_df_pivot <- pivot_longer(log2_cpm_of_DGE_df, # dataframe to be pivoted
                                  cols = "Aeration-1":"Stagnant-3", # column names to be stored as a SINGLE variable
                                  names_to = "samples", # name of that new variable (column)
                                  values_to = "expression") # name of new variable (column) storing all the values (data)

# let's look at the impact of pivoting the data (pivoted_data_plot plot)
log2_cpm_df_pivot

# Filter your data ----

# first, take a look at how many genes or transcripts have no read counts at all
table(rowSums(DGE_list$counts==0)==10)

# FALSE-9376 TRUE-5985 (true is when we have a zero in all rows aka the gene is not expressed)
# read count - counting how many times a genes has been expressed
# read = basically a sequence that gets translated into amino acids and further make proteins

# This is where the filtering starts!

# 3 aerated and 3 stagnant -> smallest group in comparison; cuttoff
keepers <- rowSums(cpm_of_DGE > 1) >= 3 

# now use base R's simple subsetting method to filter your DGE List based on the logical produced above
DGE_list_filtered <- DGE_list[keepers,]
dim(DGE_list_filtered)

log2_cpm_filtered <- cpm(DGE_list_filtered, log=TRUE)
log2_cpm_filtered_df <- as_tibble(log2_cpm_filtered, rownames = "geneID")
colnames(log2_cpm_filtered_df) <- c("geneID", sampleLabels)

# pivot this FILTERED data, just as you did earlier (pivot_longer_plot)
# pivot_longer is used to elongate the data by turning columns into rows
log2_cpm_filtered_df_pivot <- pivot_longer(log2_cpm_filtered_df, # dataframe to be pivoted
                                           cols = "Aeration-1":"Stagnant-3", 
                                           names_to = "samples", 
                                           values_to = "expression")

# Normalize your data ----
DGE_list_filtered_norm <- calcNormFactors(DGE_list_filtered, method = "TMM")

# use the 'cpm' function from EdgeR to get counts per million from your normalized data
log2_cpm_filtered_norm <- cpm(DGE_list_filtered_norm, log=TRUE)
log2_cpm_filtered_norm_df <- as_tibble(log2_cpm_filtered_norm, rownames = "geneID")
colnames(log2_cpm_filtered_norm_df) <- c("geneID", sampleLabels)

log2_cpm_filtered_norm_df_pivot <- pivot_longer(log2_cpm_filtered_norm_df, 
                                                cols = "Aeration-1":"Stagnant-3", 
                                                names_to = "samples", 
                                                values_to = "expression") 

# Identify variables of interest in study design file ----
targets
group <- targets$group
group <- factor(group)

# for the next part we are working with the normalised and filtered data (log2_cpm_filtered_norm_df_pivot)
# changing the variable name for a better overview
# using the data matrix because the clustering works only with matrices and not with data frames
norm_data <- log2_cpm_filtered_norm

# Hierarchical clustering ---------------

distance <- dist(t(norm_data), method = "euclidean")

# the dist function calculates distances only of the rows, since we are interested of clusters between groups we have to change 
# the columns to rows with the t function; safest option is to use eucledian for the method
# other agglomeration methods are "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", or "centroid"
clusters <- hclust(distance, method = "average") 

# Principal component analysis (PCA) -------------

pca_result <- prcomp(t(norm_data), scale.=F, retx=T)
pca_res_df <- as_tibble(pca_result$x)

# look at the PCA result (pca.res) that you just created
ls(pca_result)

# Prints variance summary for all principal components.
summary(pca_result) 

# $rotation shows you how much each gene influenced each principal component (called 'scores')
pca_result$rotation 

# how is each gene affected by rotation and how much does it contribute to the pc
# 'x' shows you how much each sample influenced each principal component (called 'loadings')
pca_result$x 

# sdev^2 captures these eigenvalues from the PCA result
pc_var <- pca_result$sdev^2 

# we can then use these eigenvalues to calculate the percentage variance explained by each PC
pc_percentage_variance <- round(pc_var/sum(pc_var)*100, 1)
pc_percentage_variance

# Create a PCA 'small multiples' chart ----

# this is another way to view PCA laodings to understand impact of each sample on each principal component
pca_res_df <- pca_result$x[,1:4] %>% 
  as_tibble() %>%
  add_column(sample = sampleLabels,
             group = group)

pca_pivot <- pivot_longer(pca_res_df, 
                          cols = PC1:PC4, 
                          names_to = "PC", 
                          values_to = "loadings")


# Use dplyr 'verbs' to modify our dataframe ----

# use dplyr 'mutate' function to add new columns based on existing data
mydata_df <- log2_cpm_filtered_norm_df %>% 
  mutate(aerated.AVG = (log2_cpm_filtered_norm_df$`Aeration-1` + log2_cpm_filtered_norm_df$`Aeration-2` + log2_cpm_filtered_norm_df$`Aeration-3`)/3,
         stagnant.AVG = (log2_cpm_filtered_norm_df$`Stagnant-1` + log2_cpm_filtered_norm_df$`Stagnant-2` + log2_cpm_filtered_norm_df$`Stagnant-3`)/3,
         # now make columns comparing each of the averages above that you're interested in
         LogFC = (stagnant.AVG - aerated.AVG)) %>% 
  mutate_if(is.numeric, round, 2)
mydata_df

# first, we'll use dplyr "arrange" function to sort rows based on the values in a column of interest
# then we'll display 'select' only the columns we're interested in seeing
mydata_sort <- mydata_df %>%
  dplyr::arrange(desc(LogFC)) %>% 
  dplyr::select(geneID, LogFC)
mydata_sort

# logFC --> log of fold change gives a better scale
# + fold change --> increase in expression; - fold change value --> decrease in expression
# taken to base 2 --> logFC(1) means its expressed twice as much

mydata_filter <- mydata_df %>%
  dplyr::filter(geneID=="MMP1" | geneID=="GZMB" | geneID=="IL1B" | geneID=="GNLY" | geneID=="IFNG"
                | geneID=="CCL4" | geneID=="PRF1" | geneID=="APOBEC3A" | geneID=="UNC13A" ) %>%
  dplyr::select(geneID, aerated.AVG, stagnant.AVG, LogFC) %>%
  dplyr::arrange(desc(LogFC))

# you can also filter based on any regular expression
mydata_reg_expression <- mydata_df %>%
  dplyr::filter(grepl('CXCL|IFI', geneID)) %>%
  dplyr::select(geneID, aerated.AVG, stagnant.AVG, LogFC) %>%
  dplyr::arrange(desc(geneID))

# Produce publication-quality tables using the gt package ----
gt(mydata_filter)

# Make an interactive table using the DT package ----
datatable(mydata_df[,c(1,8:10)], 
          extensions = c('KeyTable', "FixedHeader"), 
          filter = 'top',
          options = list(keys = TRUE, 
                         searchHighlight = TRUE, 
                         pageLength = 10, 
                         lengthMenu = c("10", "25", "50", "100")))

# Set up your design matrix ----

group <- factor(targets$group)

#this fits the linear model; we use 0 as start of the intercept bc in the violin plot the data was generally at the same level
design <- model.matrix(~0 + group) 
colnames(design) <- levels(group)

# Model mean-variance trend and fit linear model to data ----

# Use VOOM function from Limma package to model the mean-variance relationship
voom_DEG_list_filtered_norm <- voom(DGE_list_filtered_norm, design, plot = TRUE)

# fit a linear model to your data
mydata_linear_fit <- lmFit(voom_DEG_list_filtered_norm, design)

# the way we have first stagnant and then aerated it's like we are comparing treatment to control so at the end if we get a 
# negative value genes in stagnant group would be less expressed than in aerated
contrast_matrix <- makeContrasts(infection = stagnant - aerated,    
                                 levels=design)

# extract the linear model fit -----
mydata_linear_fit_extracted <- contrasts.fit(mydata_linear_fit, contrast_matrix)

#get bayesian stats for your linear model fit
mydata_linfit_bayes <- eBayes(mydata_linear_fit_extracted)

# TopTable to view DEGs -----
# if you change number to 10 you get the 10 highest values for logFC
mydata_DEG_topTable <- topTable(mydata_linfit_bayes, adjust ="BH", coef=1, number=40000, sort.by="logFC") 

# convert to a tibble
mydata_DEG_topTable_df <- mydata_DEG_topTable %>%
  as_tibble(rownames = "geneID")

my_data_table <- gt(mydata_DEG_topTable_df)

# TopTable (from Limma) outputs a few different stats:
# logFC, AveExpr, and P.Value should be self-explanatory
# adj.P.Val is your adjusted P value, also known as an FDR (if BH method was used for multiple testing correction)
# B statistic is the log-odds that that gene is differentially expressed. If B = 1.5, then log odds is e^1.5, where e 
# is Euler's constant (approx. 2.718).  So, the odds of differential expression is about 4.8 to 1
# t statistic is ratio of the logFC to the standard error (where the error has been moderated across all genes...
# because of Bayesian approach)

# pulling out DEG results
results <- decideTests(mydata_linfit_bayes, method="global", adjust.method="BH", p.value=0.05, lfc=1)
head(results)
summary(results)

# retrieve expression data for your DEGs ----
head(voom_DEG_list_filtered_norm$E)
colnames(voom_DEG_list_filtered_norm$E) <- sampleLabels

diffGenes <- voom_DEG_list_filtered_norm$E[results[,1] !=0,]
head(diffGenes)
dim(diffGenes) #1 genes that passed the threshold

#convert your DEGs to a dataframe using as_tibble
diffGenes_df <- as_tibble(diffGenes, rownames = "geneID")

# create interactive tables to display your DEGs ----
DEG_table <- datatable(diffGenes_df,
                       extensions = c('KeyTable', "FixedHeader"),
                       caption = 'Table 1: DEGs in Oryza sativa',
                       options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
             formatRound(columns=c(2:11), digits=2)

#write your DEGs to a file
write_tsv(diffGenes_df,"DiffGenes.txt")

# Carry out GO enrichment using gProfiler2 ----

# use topTable result to pick the top genes for carrying out a Gene Ontology (GO) enrichment analysis
mydata_DEG_topTable <- topTable(mydata_linfit_bayes, adjust ="BH", coef=1, number=50, sort.by="logFC") 

# use the 'gost' function from the gprofiler2 package to run GO enrichment analysis
gost_result <- gost(rownames(mydata_DEG_topTable), organism = "osativa", correction_method = "fdr")

# Competitive GSEA using CAMERA----

# for competitive tests the null hypothesis is that genes in the set are, at most, as often differentially 
# expressed as genes outside the set

# first let's create a few signatures to test in our enrichment analysis
signatures_real <- rownames(mydata_DEG_topTable) 
signatures_fake <- sample((rownames(voom_DEG_list_filtered_norm$E)), size = 50, replace = FALSE)
collection <- list(real = signatures_real, fake = signatures_fake)

# now test for enrichment using CAMERA
camera_result <- camera(voom_DEG_list_filtered_norm$E, collection, design, contrast_matrix[,1]) 
camera_df <- as_tibble(camera_result, rownames = "setName")
camera_df

# Self-contained GSEA using ROAST----

# remember that for self-contained the null hypothesis is that no genes in the set are differentially expressed
# mroast adjusts for multiple testing
mroast(voom_DEG_list_filtered_norm$E, collection, design, contrast=1) 

#################################################################################################################################
# WE TRIED TO ENTER THE DATA THROUGH THE CODE WE WERE PROVIDED BY THE LECTURER, BUT WE STILL COULDN'T USE THE RESULTS FOR FURTHER
# ANALYSIS!

# Quick solution for errors while importing some of the PlantGSEA gmt files
# 1. Add ".csv "extension to the downloaded file, here for rice, the file name is "Osa.DetailInfo"
# 2. Read the file
tmp = read.csv("Osa.DetailInfo.csv", header = F, sep = "\t")
# 3. make tibble
tmp = as_tibble(tmp)
# 4. remove Duplicates
tmp = tmp[!duplicated(tmp$V1), ]
# 5. write new file
write.table(tmp, "OsaUnique.csv", sep="\t",col.names = F,row.names = F)
# 6. read the file as Gmt
broadSet.Osa.Unique = getGmt("OsaUnique.csv", geneIdType=SymbolIdentifier())

#################################################################################################################################

# FROM HERE ON WE COULDN'T GET USABLE DATA, SO IT COULD BE THAT THE CODE IS PARTIALLY WRONG, AS WE WEREN'T ABLE TO COMPLETE THE
# ANALYSIS WITH THE DATA WE AQCURIED PREVIOUSLY

# camera requires collections to be presented as a list, rather than a tibble, so we must read in our signatures using the 
#'getGmt'
broadSet.C2.ALL_test <- getGmt("OSAUnique.csv", geneIdType=SymbolIdentifier())

broadSet.C2.ALL_gmt <- getGmt("OSAUnique.csv")
# before I used the file I obtained with the provided code, I used a file I generated myself and I took out one line (lipid 
# catabolic process (GO:0016042)	1	+	9.58E-01	9.58E-0) because it was a duplicate Id, so I wanted to see if it makes 
# a difference but it doesn't, I got an error that something is missing

broadSet.C2.ALL_gmt <- geneIds(broadSet.Osa.Unique)

#extract as a list
camera.res <- camera(mydata_DEG_topTable, broadSet.C2.ALL_gmt, design, contrast_matrix[,1]) 
camera.df <- as_tibble(camera.res, rownames = "setName")
camera.df


# filter based on FDR and display as interactive table
camera.df <- filter(camera.df, FDR<=0.01)

datatable(camera.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in leishmaniasis',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2,4,5), digits=2)

#as before, add a variable that maps up/down regulated pathways with phenotype
camera.df <- camera.df %>%
  mutate(phenotype = case_when(
    Direction == "Up" ~ "stagnant",
    Direction == "Down" ~ "aerated"))

#easy to filter this list based on names of signatures using 'str_detect'
#here is an example of filtering to return anything that has 'CD8' or 'CYTOTOX' in the name of the signature
camera.df.sub <- camera.df %>%
  dplyr::filter(str_detect(setName, "CD8|CYTOTOX"))

# Single sample GSEA using the GSVA package----

GSVA.res.C2CP <- gsva(voom_DEG_list_filtered_norm$E, #your data
                      broadSet.C2.ALL_gmt, #signatures
                      min.sz=3, max.sz=500, #criteria for filtering gene sets
                      mx.diff=FALSE,
                      method="gsva") #options for method are "gsva", ssgsea', "zscore" or "plage"

# Apply linear model to GSVA result
# now using Limma to find significantly enriched gene sets in the same way you did to find diffGenes
# this means you'll be using topTable, decideTests, etc
# note that you need to reference your design and contrast matrix here
fit.C2CP <- lmFit(GSVA.res.C2CP, design)
ebFit.C2CP <- eBayes(fit.C2CP)

# use topTable and decideTests functions to identify the differentially enriched gene sets
topPaths.C2CP <- topTable(ebFit.C2CP, adjust ="BH", coef=1, number=50, sort.by="logFC")
res.C2CP <- decideTests(ebFit.C2CP, method="global", adjust.method="BH", p.value=0.05, lfc=0)
# the summary of the decideTests result shows how many sets were enriched in induced and repressed genes in all sample types
summary(res.C2CP)

# pull out the gene sets that are differentially enriched between groups
diffSets.C2CP <- GSVA.res.C2CP[res.C2CP[,1] !=0,]
head(diffSets.C2CP)
dim(diffSets.C2CP)


################################################### Producing the Plots ########################################################


# producing a scatter plot of the transformed data (myTPM_stats)
my_TPM_plot <- ggplot(myTPM_stats) + 
  aes(x = SD, y = MED) +
  geom_point(shape=16, size=2) +
  geom_smooth(method=lm) +
  geom_hex(show.legend = FALSE) +
  labs(y="Median", x = "Standard deviation",
       title="Transcripts per million (TPM)",
       subtitle="unfiltered, non-normalized data",
       caption="DIYtranscriptomics - Spring 2020") +
  theme_classic() +
  theme_dark() + 
  theme_bw()
my_TPM_plot

# plotting the pivoted data
pivoted_data_plot <- ggplot(log2_cpm_df_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="unfiltered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
pivoted_data_plot

pivot_longer_plot <- ggplot(log2_cpm_filtered_df_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, non-normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
pivot_longer_plot

pivot_normalised_plot <- ggplot(log2_cpm_filtered_norm_df_pivot) +
  aes(x=samples, y=expression, fill=samples) +
  geom_violin(trim = FALSE, show.legend = FALSE) +
  stat_summary(fun = "median", 
               geom = "point", 
               shape = 95, 
               size = 10, 
               color = "black", 
               show.legend = FALSE) +
  labs(y="log2 expression", x = "sample",
       title="Log2 Counts per Million (CPM)",
       subtitle="filtered, TMM normalized",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()
pivot_normalised_plot

# producing a figure of all the pivoted data
figure_pivoted_data <- plot_grid(pivoted_data_plot, pivot_longer_plot, pivot_normalised_plot, labels = c('A', 'B', 'C'), 
                                 label_size = 12)
figure_pivoted_data

# plotting the cluster dendogram
plot(clusters, labels=sampleLabels)

# A screeplot is a standard way to view eigenvalues for each PCA
screeplot(pca_result)

# Visualize your PCA result ------------------

# We know how much each sample contributes to each PC (loadings), so let's plot
# Remember that PCA is unsupervised, so knows nothing about group assignment (aerated vs stagnant)

# But *we* know, and so we can use this knowledge to enhance the plot.  Add a 'color=group' mapping to the aes of 
# the plot above --> without this, there are only black dots

ggplot(pca_res_df) +
  aes(x = PC1, y = PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc_percentage_variance[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_percentage_variance[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  #coord_fixed() +
  theme_bw()

# Can we figure out the identity of the outlier?  We have already provided samplelabel mapping in aes, so just uncomment
# the 'geom_label()' --> from what it looks like, stagnant 1 looks like an oulier from stagnant group and aeration 3 from the 
# aeration group

ggplot(pca_res_df) +
  aes(x = PC1, y = PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc_percentage_variance[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_percentage_variance[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  #coord_fixed() +
  theme_bw()

# Uncomment 'coord_fixed()' to apply the correct aspect ratio

ggplot(pca_res_df) +
  aes(x = PC1, y = PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc_percentage_variance[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_percentage_variance[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# Uncomment 'stat_ellipse()' to see how you can circle clusters on the PCA --> couldn't calculate the ellipse bc too few points

ggplot(pca_res_df) +
  aes(x = PC1, y = PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  geom_label() +
  # stat_ellipse() +
  xlab(paste0("PC1 (",pc_percentage_variance[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_percentage_variance[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  coord_fixed() +
  theme_bw()

# How would this PCA look if you used raw counts (myCounts) instead of log2 CPM? 
# --> I don't think it works bc it's not a data frame

ggplot(myCounts) +
  aes(x = PC1, y = PC2, label=sampleLabels, color = group) +
  geom_point(size=4) +
  #geom_label() +
  #stat_ellipse() +
  xlab(paste0("PC1 (",pc_percentage_variance[1],"%",")")) + 
  ylab(paste0("PC2 (",pc_percentage_variance[2],"%",")")) +
  labs(title="PCA plot",
       caption=paste0("produced on ", Sys.time())) +
  #coord_fixed() +
  theme_bw() 

# What are the disadvantages of looking at a PCA result using such a simple XY plot?
# no idea

# plot the pivoted data frame of the PCA
ggplot(pca_pivot) +
  aes(x = sample, y = loadings, fill = group) + 
  geom_bar(stat="identity") +
  facet_wrap(~PC) +
  labs(title="PCA 'small multiples' plot",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw() +
  coord_flip()

# Make an interactive scatter plot with plotly -----

mydata_plot <- ggplot(mydata_df) +
  aes(x=aerated.AVG, y=stagnant.AVG, 
      text = paste("Symbol:", geneID)) +
  geom_point(shape=16, size=1) +
  ggtitle("Stagnant vs. Aerated") +
  theme_bw()

mydata_interactive <- ggplotly(mydata_plot)

# Volcano Plots ----

# Volcano plot of my Top table hits (values at top are the most significant)
vplot_of_mydata <- ggplot(mydata_DEG_topTable_df) +
  aes(y=-log10(adj.P.Val), x=logFC, text = paste("Symbol:", geneID)) +
  geom_point(size=2) +
  geom_hline(yintercept = -log10(0.01), linetype="longdash", colour="grey", linewidth=1) +
  geom_vline(xintercept = 1, linetype="longdash", colour="#BE684D", linewidth=1) +
  geom_vline(xintercept = -1, linetype="longdash", colour="#2C467A", linewidth=1) +
  labs(title="Volcano plot",
       subtitle = "Oryza sativa",
       caption=paste0("produced on ", Sys.time())) +
  theme_bw()

ggplotly(vplot_of_mydata)

vennDiagram(results, include="both")

# Create a heatmap of differentially expressed genes ----

heatmaply(diffGenes_df[1:50,2:7], #added the 1:50 to make it easier to visualise (just the first 50 genes)
          xlab = "Samples", ylab = "DEGs", 
          main = "DEGs in Oryza sativa",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = FALSE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffGenes_df)[2:7],
          labRow = diffGenes_df$geneID[1:50],
          heatmap_layers = theme(axis.line=element_blank())
)


# produce an interactive manhattan plot of enriched GO terms
# set interactive = FALSE to get plot for publications
mydata_gostplot <- gostplot(gost_result, interactive = F, capped = F) 

# produce a publication quality static manhattan plot with specific GO terms highlighted
# rerun the above gostplot function with 'interactive=F' and save to an object 'mygostplot'
publish_gostplot(
  mydata_gostplot,
  # highlight_terms = c("GO:0019825"),
  filename = NULL,
  width = NA,
  height = NA)

gostplot(gost_result, interactive = T, capped = F)

##############################################################################################################################
# THE FOLLOWING CODE WAS NOT ADJUSTED TO OUR DATA, AS THE COMPETITIVE ANALYSIS DIDN'T WORK AND THERE WAS NO DATA TO USE FURTHER

# graph camera results as bubble chart 
ggplot(camera.df[1:25,], aes(x=phenotype, y=setName)) + 
  geom_point(aes(size=NGenes, color = bluered(n), alpha=-log10(FDR))) +
  theme_bw()

# Create a heatmap of differentially expressed genes ----

heatmaply(diffSets.C2CP, 
          #dendrogram = "row",
          xlab = "Samples", ylab = "KEGG pathways", 
          main = "Responsive KEGG pathways in cutaneous leishmaniasis",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.0000001,
          titleX = T,
          titleY = T,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          label_names = c("Gene", "Sample:", "Value"),
          fontsize_row = 5, fontsize_col = 5,
          labCol = colnames(diffSets.C2CP),
          labRow = rownames(diffSets.C2CP),
          heatmap_layers = theme(axis.line=element_blank())
)