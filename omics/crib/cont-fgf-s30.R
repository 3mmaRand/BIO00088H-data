
# OMICS 1 HELLO DATA ------------------------------------------------------

# Load tidyverse 
library(tidyverse)


# Import ------------------------------------------------------------------

# 🐸 import the s30 data
s30 <- read_csv("omics/week-3/data-raw/xlaevis_counts_s30.csv")


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
          file = "omics/week-3/data-processed/s30_summary_samp.csv")

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
          file = "omics/week-3/data-processed/s30_summary_gene.csv")


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
          file = "omics/week-3/data-processed/s30_summary_gene_filtered.csv")


# Use the list of `xenbase_gene_id` in the filtered summary to filter
# the original dataset:
s30_filtered <- s30 |> 
  filter(xenbase_gene_id %in%  s30_summary_gene_filtered$xenbase_gene_id)

# Write the filtered original to file:
write_csv(s30_filtered, 
          file = "omics/week-3/data-processed/s30_filtered.csv")



# OMICS 2 STATISTICAL ANALYSIS --------------------------------------------
library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

library(DESeq2)

# import the filtered s30
s30_filtered <- read_csv("omics/week-3/data-processed/s30_filtered.csv")

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
write_csv(s30_fgf_only, "omics/week-4/results/s30_fgf_only.csv")


## Import metadata that maps the sample names to treatments
meta <- read_table("omics/week-5/meta/frog_meta_data.txt")
row.names(meta) <- meta$sample_id

# We only need the s30
meta_s30 <- meta |>
  dplyr::filter(stage == "stage_30")



# 1. check the names match: rownames of meta and column names of counts
names(s30_filtered)
row.names(meta_s30)

# verfiy with equality
names(s30_filtered[2:7]) == row.names(meta_s30)



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
                              colData = meta_s30,
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


# Notice there is a relationship between the total number of 
# counts and the size factors which we can see with a quick scatter plot
plot(colSums(s30_count_mat), sizeFactors(dds))

# Then to get the normalized counts as a matrix, we can use:
normalised_counts <- counts(dds, normalized = TRUE)


# Save this normalized data matrix to file for later use:
data.frame(normalised_counts,
           xenbase_gene_id = row.names(normalised_counts)) |>
  write_csv(file = "omics/week-4/results/s30_normalised_counts.csv")

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
  write_csv(file = "omics/week-4/results/s30_results.csv")



# OMICS 3 VISUALISE AND INTERPRET  -----------------------------------------

library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")

# Import the normalised counts
s30_count_norm <- read_csv("omics/week-4/results/s30_normalised_counts.csv")

s30_count_norm

# DE results
s30_results <- read_csv("omics/week-4/results/s30_results.csv")
s30_results


# merge the results with the normalised counts
s30_results <- s30_count_norm |>
  left_join(s30_results, by = "xenbase_gene_id")


# LINKING INFORMATION -----------------------------------------------------

# This information comes from the xenbase information pages

# xenbase Gene Product Information [readme]gzipped gpi (tab separated)
# File contains gene product information for species specific Xenbase genes. 
# This file is in the Gene Product information 2.1 format and is provided with 
# gzip compression.
# https://download.xenbase.org/xenbase/GenePageReports/xenbase.gpi.gz
# gunzip xenbase.gpi.gz
# less xenbase.gpi
# q
# This file contains gene product information for both 
# Xenopus tropicalis (taxon:8364) and Xenopus laevis (taxon:8355)

library(readxl)
gene_info <- read_excel("omics/week-5/meta/xenbase_info.xlsx") 

# join the gene info with the results
s30_results <- s30_results |>
  left_join(gene_info, by = "xenbase_gene_id")


## Import metadata that maps the sample names to treatments
meta <- read_table("omics/week-4/meta/frog_meta_data.txt")
row.names(meta) <- meta$sample_id

# We only need the s30
meta_s30 <- meta |>
  dplyr::filter(stage == "stage_30")



# log2 transformed normalised counts needed for data viz

# log2 transform the counts in s30_count_norm
# add a tiny amount to avoid log(0)
s30_results <- s30_results |>
  mutate(across(starts_with("s30"), 
                \(x) log2(x + 0.001),
                .names = "log2_{.col}"))

# https://posit-conf-2023.github.io/programming-r/03-iteration-01.html#/title-slide

# We now have datatframe with all the info we need, normalised counts, 
# log2 normalised counts, statistical comparisons with fold changes and p values, 
# information about the gene
# other than just the id


s30_results_sig0.01 <- s30_results |> 
  filter(padj <= 0.01)
# 59 genes

# write to csv file
write_csv(s30_results_sig0.01, 
          file = "omics/week-5/results/s30_results_sig0.01.csv")

s30_results_sig0.05 <- s30_results |> 
  filter(padj <= 0.05)

# 117 genes

# write to csv file
write_csv(s30_results_sig0.05, 
          file = "omics/week-5/results/s30_results_sig0.05.csv")



# PCA AND CLUSTERING ------------------------------------------------------


# we do this on the log2 transformed normalised counts or the regularized the
# log transformed counts

# transpose the data
# we are reducing he number of dimensions from 10136
s30_log2_trans <- s30_results |> 
  select(starts_with("log2_")) |>
  t() |> 
  data.frame()

colnames(s30_log2_trans) <- s30_results$xenbase_gene_id

# just for indep study before
# s30_log2_trans$sample <- row.names(s30_log2_trans)
# a <- s30_log2_trans |> ggplot(aes(x = `XB-GENE-1000007`, 
#                              y = `XB-GENE-1000023`)) +
#   geom_point() +
#   geom_text(aes(label = sample), 
#             vjust = -1, size = 3) +
#   scale_x_continuous(expand = c(0.05,0.05)) +
#   scale_y_continuous(expand = c(0.05,0.05)) +
#   theme_classic()
# 
# 
# b <- s30_log2_trans |> ggplot(aes(x = `XB-GENE-1000062`, 
#                                   y = `XB-GENE-1000072`)) +
#   geom_point() +
#   geom_text(aes(label = sample), 
#             vjust = -1, size = 3) +
#   scale_x_continuous(expand = c(0.05,0.05)) +
#   scale_y_continuous(expand = c(0.05,0.05)) +
#   theme_classic()
# 
# c <- s30_log2_trans |> ggplot(aes(x = `XB-GENE-1000113`, 
#                                   y = `XB-GENE-1000132`)) +
#   geom_point() +
#   geom_text(aes(label = sample), 
#             vjust = -1, size = 3) +
#   scale_x_continuous(expand = c(0.05,0.05)) +
#   scale_y_continuous(expand = c(0.05,0.05)) +
#   theme_classic()
# 
# d <- s30_log2_trans |> ggplot(aes(x = `XB-GENE-1000149`, 
#                                   y = `XB-GENE-1000251`)) +
#   geom_point() +
#   geom_text(aes(label = sample), 
#             vjust = -1, size = 3) +
#   scale_x_continuous(expand = c(0.05,0.05)) +
#   scale_y_continuous(expand = c(0.05,0.05)) +
#   theme_classic()
# 
# library(patchwork)
# fig <- (a + b) / (c + d)
# 
# ggsave("omics/week-5/images/why_pca.png", 
#        plot = fig,
#        width = 6, height = 6)

# perform PCA using standard functions
pca <- s30_log2_trans |>
  prcomp(scale. = TRUE) 

summary(pca)
# Importance of components:
#                         PC1     PC2     PC3     PC4     PC5       PC6
# Standard deviation     57.1519 50.3285 43.8382 35.4644 33.9440 2.403e-13
# Proportion of Variance  0.3224  0.2500  0.1897  0.1241  0.1137 0.000e+00
# Cumulative Proportion   0.3224  0.5724  0.7621  0.8863  1.0000 1.000e+00

# remove log2 from  row.names(s30_log2_trans)
# to labelthe pca results
sample_id <- row.names(s30_log2_trans) |> str_remove("log2_")

pca_labelled <- data.frame(pca$x,
                           sample_id)

# merge with metadata
# so we can label points by treatment and sib group
pca_labelled <- pca_labelled |> 
  left_join(meta_s30, 
            by = "sample_id")

pca <- pca_labelled |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = sibling_rep,
             shape = treatment)) +
  geom_point(size = 3) +
  scale_colour_viridis_d(end = 0.95, begin = 0.15,
                         name = "Sibling pair",
                         labels = c("A", ".5", ".6")) +
  scale_shape_manual(values = c(21, 19),
                     name = NULL,
                     labels = c("Control", "FGF-Treated")) +
  theme_classic()

# There is a bit of separation between treatments on PCA2
# not that it isn't easy to draw strong conclusions on the basis of 3 points

ggsave("omics/week-5/figures/frog-s30-pca.png",
       plot = pca,
       height = 3, 
       width = 4,
       units = "in",
       device = "png")




# Heatmap
# only should do on sig genes.
# but use the log 2 normalised values

mat <- s30_results_sig0.01 |> 
  select(starts_with("log2_")) |>
  as.matrix()

rownames(mat) <- s30_results_sig0.01$xenbase_gene_symbol

n_treatment_clusters <- 2
n_gene_clusters <- 2

library(heatmaply)

heatmaply(mat, 
          scale = "row",
          hide_colorbar = TRUE,
          k_col = n_treatment_clusters,
          k_row = n_gene_clusters,
          label_names = c("Gene", "Sample", "Expression (normalised, log2)"),
          fontsize_row = 7, fontsize_col = 10,
          labCol = str_remove(colnames(mat), pattern = "log2_"),
          labRow = rownames(mat),
          heatmap_layers = theme(axis.line = element_blank()))


# volcano plot
# s30_results

# colour the points if padj < 0.05
# and log2FoldChange > 1
library(ggrepel)


s30_results <- s30_results |> 
  mutate(log10_padj = -log10(padj),
         sig = padj < 0.05,
         bigfc = abs(log2FoldChange) >= 2) 

vol <- s30_results |> 
  ggplot(aes(x = log2FoldChange, 
             y = log10_padj, 
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed") +
  geom_vline(xintercept = 2, 
             linetype = "dashed") +
  geom_vline(xintercept = -2, 
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("gray", 
                                 "pink",
                                 "gray30",
                                 "deeppink")) +
  geom_text_repel(data = subset(s30_results, 
                                bigfc & sig),
                  aes(label = xenbase_gene_symbol),
                  size = 3,
                  max.overlaps = 50) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("omics/week-5/figures/frog-s30-volcano.png",
       plot = vol,
       height = 4.5, 
       width = 4.5,
       units = "in",
       device = "png")

