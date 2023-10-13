
# OMICS 1 HELLO DATA ------------------------------------------------------

# Load tidyverse 
library(tidyverse)


# Import ------------------------------------------------------------------

# 🐸 import the s30 data
s30 <- read_csv("data-raw/xlaevis_s30_filtered.csv")


# Explore -----------------------------------------------------------------


# Distribution of values across the whole dataset -------------------------

# raw values
s30 |>
  pivot_longer(cols = -xenbase_gene_id,
               names_to = "sample",
               values_to = "count") |>
  ggplot(aes(x = count)) +
  geom_histogram()

# This data is very skewed - there are so many low values that we can't
# see the tiny bars for the higher values. Logging the counts is a way to
# make the distribution more visible.

# Repeat the plot on log of the counts.
s30 |>
  pivot_longer(cols = -xenbase_gene_id,
               names_to = "sample",
               values_to = "count") |>
  ggplot(aes(x = log10(count))) +
  geom_histogram()

# peak at zero suggests quite a few counts of 1. We would expect we
# would expect the distribution of counts to be roughly log normal 
# because this is expression of *all* the genes in the genome
# small peak near the low end suggests that these lower counts might be anomalies.
# The excess number of low counts indicates we might want to create a cut
# off for quality control. The removal of low counts is a common
# processing step in 'omic data. We will revisit this after we have
# considered the distribution of counts across samples and genes.


#  Distribution of values across the sample/cells -------------------------
# Get a quick overview of the columns:
summary(s30)

# Notice that: - the minimum count is 0 and the maximums are very high in
# all the columns - the medians are quite a lot lower than the means so
# the data are skewed (hump to the left, tail to the right) - there must
# be quite a lot of zeros - the columns are roughly similar and it doesn't
# look like there is an anomalous replicate.

# Summarise all the samples:
s30_summary_samp <- s30 |>
  pivot_longer(cols = -xenbase_gene_id,
               names_to = "sample",
               values_to = "count") |>
  group_by(sample) |>
  summarise(min = min(count),
            lowerq = quantile(count, 0.25),
            mean = mean(count),
            median = median(count),
            upperq = quantile(count, 0.75),
            max = max(count),
            n_zero = sum(count == 0))
# The mean count ranges from 260 to 426.

# Write `s30_summary_samp` to a file called "s30_summary_samp.csv":
write_csv(s30_summary_samp, 
          file = "data-processed/s30_summary_samp.csv")

# Pivot the data and plot logged counts
s30 |>
  pivot_longer(cols = -xenbase_gene_id,
               names_to = "sample",
               values_to = "count") |>
  ggplot(aes(log10(count))) +
  geom_density() +
  facet_wrap(. ~ sample, nrow = 3)

# The key information to take from these plots is:
# # the distributions are roughly similar in width, height, location and
#     overall shape so it doesn't look as though we have any suspect
# samples
# the peak at zero suggests quite a few counts of 1.
# since we would expect the distribution of counts in each sample to
# be roughly log normal so that the small rise near the low end, even
# before the peak at zero, suggests that these lower counts might be
# anomalies.
# 
# The excess number of low counts indicates we might want to create a cut
# off for quality control. The removal of low counts is a common
# processing step in 'omic data. We will revisit this after we have
# considered the distribution of counts across genes (averaged over the
# samples).


# Distribution of values across the genes ---------------------------------

# Summarise the counts for each genes:
s30_summary_gene <- s30 |>
  pivot_longer(cols = -xenbase_gene_id,
               names_to = "sample",
               values_to = "count") |>
  group_by(xenbase_gene_id) |>
  summarise(min = min(count),
            lowerq = quantile(count, 0.25),
            sd = sd(count),
            mean = mean(count),
            median = median(count),
            upperq = quantile(count, 0.75),
            max = max(count),
            total = sum(count),
            n_zero = sum(count == 0))

# Notice that we have:
# 
# a lot of genes with counts of zero in every sample
# a lot of genes with zero counts in several of the samples
# some very very low counts.
# 
# These should be filtered out because they are unreliable - or, at the
# least, uninformative. The goal of our downstream analysis will be to see
# if there is a significant difference in gene expression between the
# control and FGF-treated sibling. Since we have only three replicates in
# each group, having one or two unreliable, missing or zero values, makes
# such a determination impossible for a particular gene. We will use the
# total counts and the number of samples with non-zero values to filter
# our genes later.
# 
# As we have a lot of genes, it is again helpful to plot the mean counts
# with pointrange to get an overview. We will plot the log of the counts -
# we saw earlier that logging made it easier to understand the
# distribution of counts over such a wide range. We will also order the
# genes from lowest to highest mean count.

# Plot the logged mean counts for each gene in order of size 
s30_summary_gene |> 
  ggplot(aes(x = reorder(xenbase_gene_id, mean), y = log10(mean))) +
  geom_pointrange(aes(ymin = log10(mean - sd), 
                      ymax = log10(mean + sd )),
                  size = 0.1)

# You can see we also have quite a few genes with means less than 1 (log
# below zero). Note that the variability between genes (average counts
# between 0 and 102586) is far greater than between samples (average
# counts from 260 to 426) which is exactly what we would expect to see.

# Write `s30_summary_gene` to a file called "s30_summary_gene.csv":
write_csv(s30_summary_gene, 
          file = "data-processed/s30_summary_gene.csv")


# Filtering for QC --------------------------------------------------------

# Our samples look to be similarly well sequenced. There are no samples we
# should remove. However, some genes are not express or the expression
# values are so low in for a gene that they are uninformative. We will
# filter the `s30_summary_gene` dataframe to obtain a list of
# `xenbase_gene_id` we can use to filter `s30`.

# My suggestion is to include only the genes with counts in at least 3
# samples[^3] and those with total counts above 20.

# Filter the summary by gene dataframe:
s30_summary_gene_filtered <- s30_summary_gene |> 
  filter(total > 20) |> 
  filter(n_zero < 4)
# this leaves 10136 genes


# Write the filtered summary by gene to file:
write_csv(s30_summary_gene_filtered, 
          file = "data-processed/s30_summary_gene_filtered.csv")


# Use the list of `xenbase_gene_id` in the filtered summary to filter
# the original dataset:
s30_filtered <- s30 |> 
  filter(xenbase_gene_id %in%  s30_summary_gene_filtered$xenbase_gene_id)

# Write the filtered original to file:
write_csv(s30_filtered, 
          file = "data-processed/s30_filtered.csv")



# OMICS 2 STATISTICAL ANALYSIS --------------------------------------------
library(tidyverse)
library(conflicted)
library(DESeq2)

# import the filtered s30
s30_filtered <- read_csv("data-processed/s30_filtered.csv")

# number of genes and samples. remember we filtered out any genes 
# with 4, 5 or 6 zeros and those where the total count was less than 20


# first find genes that are expressed only in one group. that is those
# only zeros in one group but values in all of the others.

# Find the genes that are expressed only in the FGF-treated group:
  
s30_fgf_only <- s30_filtered |> 
  filter(S30_C_5 == 0, 
         S30_C_6 == 0, 
         S30_C_A == 0, 
         S30_F_5 > 0, 
         S30_F_6 > 0, 
    S30_F_A > 0)
# there are 26 genes expressed in the FGF-treated group but not in the control group

# genes that are expressed only in the control group.

s30_con_only <- s30_filtered |> 
  filter(S30_C_5 > 0, 
         S30_C_6 > 0, 
         S30_C_A > 0, 
         S30_F_5 == 0, 
         S30_F_6 == 0, 
         S30_F_A == 0)
# there are no genes expressed in the control group but not in the FGF-treated group

# Write to file 
write_csv(s30_fgf_only, "results/s30_fgf_only.csv")


## Import metadata that maps the sample names to treatments
meta <- read_table("meta/frog_meta_data.txt")
row.names(meta) <- meta$sample_id

# We only need the s30
meta_S30 <- meta |>
  dplyr::filter(stage == "stage_30")



# 1. check the names match: rownames of meta and column names of counts
names(s30_filtered)
row.names(meta_S30)

# verfiy with equality
names(s30_filtered[2:7]) == row.names(meta_S30)



# 2. Create a DESeq2 object
# bioconductor has custom data types
# to make a DESeqDataSet  object we need 
#    a) count matrix, 
#    b) metadata and
#    c) a design matrix

# The count matrix must contain only the counts
# the genes are not inside the matrix, but are the row names
s30_count_mat <- s30_filtered |>
  select(-xenbase_gene_id) |>
  as.matrix()
# add the genes as row names
row.names(s30_count_mat) <- s30_filtered$xenbase_gene_id

View(s30_count_mat)
knitr::kable(head(s30_count_mat))


# creates the DESeqDataSet (dds) object
dds <- DESeqDataSetFromMatrix(s30_count_mat,
                              colData = meta_S30,
                              design = ~ treatment + sibling_rep)
# design matrix is a formula that describes the experimental design using the 
# column names of the metadata dataframe


# This is fine: 
# converting counts to integer mode
# Warning message:

#   some variables in design formula are characters, converting to factors



# To see the counts in the DESeqDataSet object, we can use the counts() function
counts(dds) |> View()
# This should be what is in s30_count_mat

# The normalisation is done automatically with the DE
# However, we need the normalised counts for data viz.


#  The normalised counts can be generated using:
#  estimateSizeFactors() for estimating the factors for normalisation.
dds <- estimateSizeFactors(dds)
# By assigning the results back to the dds object, we are filling in the
# slots of the DESeqDataSet object with the appropriate information.

# We can take a look at the normalization factors of each sample using:
sizeFactors(dds)
# S30_C_5   S30_C_6   S30_C_A   S30_F_5   S30_F_6   S30_F_A 
# 0.8812200 0.9454600 1.2989886 1.0881870 1.0518961 0.8322894 

# Notice there is a relationship between the total number of 
# counts and the size factors which we can see with a quick scatter plot
plot(colSums(s30_count_mat), sizeFactors(dds))

# Then to get the normalized counts as a matrix, we can use:
normalised_counts <- counts(dds, normalized = TRUE)


# Save this normalized data matrix to file for later use:
data.frame(normalised_counts,
           xenbase_gene_id = row.names(normalised_counts)) |>
  write_csv(file = "results/S30_normalised_counts.csv")

# These normalized counts will be useful for downstream visualization
# of results, but cannot be used as input to DESeq2 or any other tools that
# perform differential expression analysis that use the negative binomial model.


# DE WITH DESeq2 ----------------------------------------------------------

# run the actual differential expression analysis,
# we use a single call to the function DESeq().
# Run analysis
dds <- DESeq(dds)

# check dispersion estimates, checks the assumption of the DESeq2 model
plotDispEsts(dds)
# look fine


## Define contrasts for FGF overexpression
contrast_fgf <- c("treatment", "FGF", "control")
results_fgf <- results(dds,
                       contrast = contrast_fgf,
                       alpha = 0.01)
results_fgf |>
  data.frame() |> head()


data.frame(results_fgf,
           xenbase_gene_id = row.names(results_fgf)) |>
  write_csv(file = "results/S30_results.csv")












# OMICS 3 VISUALISE AND INTERPRET and chek the assumptions?? -----------------------------------------



# PCA AND CLUSTERING -----------------------

# we do this on the log2 transformed normalised counts or the regularized the
# log transformed counts
# rlog is a method to bias from low count genes. 
# https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/03_DGE_QC_analysis.html gives a good explanation of regularized the log transform (rlog)


# The rlog transformation of the normalized counts is only necessary
# for these visualization methods during this quality assessment.
# They are not used for DE because DESeq2 takes care of that

### Transform counts for data visualization
rld <- rlog(dds, blind = TRUE)
# transformed counts are accessed with
# assay(rld)
rlogged <- assay(rld) |> t()
# These techniques want cells (samples in rows) and genes in columns.


# perform PCA using standard functions
pca <- rlogged |>
  prcomp(scale. = TRUE) 

summary(pca)
# PC1 proportion 0.3857, PC2 proportion 0.2112

pca_labelled <- data.frame(pca$x,
                           sample_id = row.names(rlogged))


# merge with metadata
pca_labelled <- pca_labelled |> 
  left_join(meta_S30, 
            by = "sample_id")

pca_labelled |> ggplot(aes(x = PC1, y = PC2, 
                           colour = treatment,
                           shape = sibling_rep)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(end = 0.95) +
  theme_classic()

# they don't cluster that well, by treatment or sib group


# ########## TO COPY TO OTHER COMPARISONS FROM HERE
# # Hierarchical Clustering
# #
# #  pheatmap need a matrix or dataframe so use SummarizedExperiment::assay()
# #  which is loaded with DESeq2
# rld_mat <- assay(rld)
# 
# # Compute pairwise correlation values
# rld_cor <- cor(rld_mat)    ## cor() is a base R function
# 
# meta_S30 <- as.data.frame(meta_S30)
# row.names(meta_S30) <-  meta_S30$sample_id
# # Plot heatmap using the correlation matrix and the metadata object
# pheatmap(rld_cor,
#          annotation = meta_S30[3:4])



## Principle Components Analysis (PCA) 
  
# explore data further - at the sample and gene level check reps cluster 
# together do on the normalised counts


# venn diagram for frogs, use original dataset

# look up information about the genes

# volcano plots

# MAYBE WEEK 3
# plotMA(results_fgf)
# 
# resultsNames(dds)
# # "Intercept", "treatment_FGF_vs_control" "sibling_rep_five_vs_A"    "sibling_rep_six_vs_A"
# 
# 
# 
# # results_fgf_shrunken <- lfcShrink(dds,
# #                                   coef = "treatment_FGF_vs_control",
# #                                   type = "apeglm")
# # 
# # plotMA(results_fgf_shrunken)