# Load pckages ------------------------------------------------------------ 
library(tidyverse)
library(scran)
library(conflicted)

# OMICS 1 HELLO DATA ------------------------------------------------------

# HSPC --------------------------------------------------------------------


# Import ------------------------------------------------------------------

# 🐭 import the hspc data
hspc <- read_csv("omics/week-3/data-raw/surfaceome_hspc.csv")

# Explore -----------------------------------------------------------------


# Distribution of values across the whole dataset -------------------------

hspc |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |> 
  ggplot(aes(x = expr)) +
  geom_histogram()

# This is a very striking distribution. Is it what we are expecting?
# Again,the excess number of low values is almost certainly anomalous.
# They will be inaccurate measure and we will want to exclude expression
# values below (about) 1. We will revisit this after we have considered
# the distribution of expression across cells and genes.
# 
# bimodal because this is a subset of the genome and the nature of the 
# subsetting has had an influence here.
# These are a subset of cell surface proteins that show a significant
# difference between at least two of twelve cell subtypes. That is, all of
# these genes are either high or low.

#  Distribution of values across the sample/cells -------------------------

# Get a quick overview the first 20 columns:
summary(hspc[1:20])

# data are logged (to base 2) normalised counts, not raw counts as 
# a minimum value of 0 appears in all 20 columns - perhaps that is
# true across the whole dataset (or at least common)
# at least some of the medians are zeros so there must be quite a lot
# of zeros
# the few columns we can see are roughly similar
# it would not be very practical to plot the distributions of values
# in cell cell using `facet_wrap()`.

# Summarise all the cells:
hspc_summary_samp <- hspc |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(cell) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            mean = mean(expr),
            median = median(expr),
            sd = sd(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            n_zero = sum(expr == 0))

# All cells have quite a few zeros and the lower quartile is 0 for all
# cells, *i.e.*, every cell has many genes with zero expression.

# Create an ordered pointrange plot.
hspc_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

# the average expression in cells is similar for all cells. This is
#     good to know - if some cells had much lower expression perhaps there
#     is something wrong with them, or their sequencing, and they should
#     be excluded.
# the distributions are roughly similar in width too

# Write `hspc_summary_samp` to a file called "hspc_summary_samp.csv":
write_csv(hspc_summary_samp, 
          file = "omics/week-3/data-processed/hspc_summary_samp.csv")


# Distribution of values across the genes ---------------------------------

# Summarise the expression for each genes:
hspc_summary_gene <- hspc |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(ensembl_gene_id) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            sd = sd(expr),
            mean = mean(expr),
            median = median(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            total = sum(expr),
            n_zero = sum(expr == 0))

# no genes with 0 in every cell
# very few genes (9) with no zeros at all
# quite a few genes with zero in many cells but this matters less than
#     zeros in the frog samples because we had just 6 samples and we have
#     701 cells.

# Plot the logged mean counts for each gene in order of size using
hspc_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd),
                  size = 0.1)
# variability between genes (average expression between 0.02 and and 10.03) 
# is far greater than between cells (average expression from1.46 to 3.18)
# which is expected.

# Write `hspc_summary_gene` to a file called "hspc_summary_gene.csv":
write_csv(hspc_summary_gene, 
          file = "omics/week-3/data-processed/hspc_summary_gene.csv")



# Filtering for QC --------------------------------------------------------

# We will take a different approach to filtering the single cell data. For
# the Frog samples we are examining the control and the FGF treated
# samples. This means have a low number of counts overall means the gene
# is not really expressed (detected) in any condition, and filtering out
# those genes is removing things that definitely are not interesting. For
# the mice, we have examined only one cell type but will be making
# comparisons between cells types. It may be that low expression of a gene
# in this cell type tells us something if that gene is highly expressed in
# another cell type. Instead, we will make statistical comparisons between
# the cell types and then filter based on overall expression, the
# difference in expression between cell types and whether that difference
# is significant.
# 
# The number of "replicates" is also important. When you have only three
# in each group it is not possible to make statistical comparisons when
# several replicates are zero. This is less of an issue with single cell
# data.


# PROG --------------------------------------------------------------------

# Import ------------------------------------------------------------------

# 🐭 import the prog data
prog <- read_csv("omics/week-3/data-raw/surfaceome_prog.csv")

# Explore -----------------------------------------------------------------


# Distribution of values across the whole dataset -------------------------

prog |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |> 
  ggplot(aes(x = expr)) +
  geom_histogram()

# bimodal like hspc

#  Distribution of values across the sample/cells -------------------------

# Get a quick overview the first 20 columns:
summary(prog[1:20])

# like hspc

# Summarise all the cells:
prog_summary_samp <- prog |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(cell) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            mean = mean(expr),
            median = median(expr),
            sd = sd(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            n_zero = sum(expr == 0))

# All cells have quite a few zeros and the lower quartile is 0 for all
# cells, *i.e.*, every cell has many genes with zero expression.

# Create an ordered pointrange plot.
prog_summary_samp |> 
  ggplot(aes(x = reorder(cell, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd ),
                  size = 0.1)

# the average expression in cells is similar for all cells. like hspc

# Write `prog_summary_samp` to a file called "prog_summary_samp.csv":
write_csv(prog_summary_samp, 
          file = "omics/week-3/data-processed/prog_summary_samp.csv")


# Distribution of values across the genes ---------------------------------

# Summarise the expression for each genes:
prog_summary_gene <- prog |>
  pivot_longer(cols = -ensembl_gene_id,
               names_to = "cell",
               values_to = "expr") |>
  group_by(ensembl_gene_id) |>
  summarise(min = min(expr),
            lowerq = quantile(expr, 0.25),
            sd = sd(expr),
            mean = mean(expr),
            median = median(expr),
            upperq = quantile(expr, 0.75),
            max = max(expr),
            total = sum(expr),
            n_zero = sum(expr == 0))

# no genes with 0 in every cell
# very few genes (12) with no zeros at all
# quite a few genes with zero in many cells but this matters less than
#     zeros in the frog samples because we had just 6 samples and we have
#     799 cells. Similar to hspc

# Plot the logged mean counts for each gene in order of size using
prog_summary_gene |> 
  ggplot(aes(x = reorder(ensembl_gene_id, mean), y = mean)) +
  geom_pointrange(aes(ymin = mean - sd, 
                      ymax = mean + sd),
                  size = 0.1)
# variability between genes (average expression between 0.02 and and 10.31) 
# is far greater than between cells (average expression from 1.32 to 3.13)
# which is expected.

# Write prog_summary_gene to a file called "prog_summary_gene.csv":
write_csv(prog_summary_gene, 
          file = "omics/week-3/data-processed/prog_summary_gene.csv")




# OMICS 2 STATISTICAL ANALYSIS --------------------------------------------

# We will carry out several steps
# 
# 1.  the data should be imported already
# 2.  Combine the two datasets ready for analysis
# 3.  Filter the data to remove genes that are not expressed in any cell
# 4.  Find the genes that are expressed in only one cell type 
# (the prog or the hspc) 
# 5.  Do differential expression analysis on the genes using the **`scran`** package. 
# This needs to be done on the logged normalised counts.

# 2. Combine the two datasets ready for analysis ---------------------------

# We need to combine the two datasets of 701 and 798 cells into one dataset 
# of 1499 cells, i.e., 1499 columns. The number of rows is the number of genes,
# 280. Before combining, we must make sure genes in the same order in 
# both dataframes or we would be comparing the expression of one gene 
# in one cell type to the expression of a different gene in the other cell type!
  
# Check the gene ids are in the same order in both dataframes:
identical(prog$ensembl_gene_id, hspc$ensembl_gene_id)
# the names are the same and in the same order

# **`scran`** can use a matrix or a dataframe of counts but theses must be 
# log normalised counts. If using a dataframe, the columns must only 
# contain the expression values (not the gene ids).

# Combine the two dataframes (minus the gene ids) into one dataframe 
# called `prog_hspc`:
prog_hspc <- bind_cols(prog[-1], hspc[-1])

# Now add the gene ids as the row names:
row.names(prog_hspc) <- prog$ensembl_gene_id

# Filter to remove unexpressed genes --------------------------------------
# In this dataset, we will not see and genes that are not expressed in any of the 
# cells because we are using a specific subset of the transcriptome that was 
# deliberately selected. However, we will go through how to do this because
# it is an important step in most analyses. 
# 
# For the 🐸 frog data you should remember that we were able to filter 
# out our unexpressed genes in [Omics 1](../week-3/workshop.html) because 
# we were examining both groups to be compared. In that workshop, 
# we discussed that we could not filter out unexpressed genes in 
# the 🐭 mouse data because we only had one cell types at that time. 
# During the Consolidate Independent Study you examined the hspc cells. 

# Where the sum of all the values in the rows is zero, all the entries must be 
# zero. We can use this to find the filter the genes that are not expressed 
# in any of the cells. To do row wise aggregates such as the sum across rows 
# we can use the `rowwise()` function. `c_across()` allows us to use the 
# colon notation `Prog_001:HSPC_852` in `sum()` rather than having to list all 
# the column names: `sum(Prog_001, Prog_002, Prog_002, Prog_004,.....)` 

# Find the genes that are 0 in every column of the prog_hspc dataframe:

prog_hspc |> 
  rowwise() |> 
  dplyr::filter(sum(c_across(Prog_001:HSPC_852)) == 0)
# There are no genes that are completely unexpressed in this set of 280 genes

# We might also examine the genes which are least expressed.
# Find ten least expressed genes:
rowSums(prog_hspc) |> sort() |> head(10)


# When you consider that there are 1499 cells, a values of 30 are low even  
# considering these are already logged and normalised (ie., the range of  
# values is less that it would be for raw counts) 
  
# Find the genes that are expressed in only one cell type -----------------
# To find the genes that are expressed in only one cell type, 
# we can use the same approach as above but only sum the columns for one cell type. 
# 
# Find the genes that are 0 in every column for the HSPC cells:
  
prog_hspc |> 
  rowwise() |> 
  dplyr::filter(sum(c_across(HSPC_001:HSPC_852)) == 0)

# Note that if we knew there were some rows that were all zero across both 
# cell types, we would need to add 
# |> dplyr::filter(sum(c_across(Prog_001:Prog_852)) != 0)

# Genes that are 0 in every column for the Prog cells:
prog_hspc |> 
  rowwise() |> 
  dplyr::filter(sum(c_across(Prog_001:Prog_852)) == 0)
# there are no genes that are expressed in only one cell type 
  
# Differential expression analysis ----------------------------------------

# Like **`DESeq2`**, **`scran`** uses a statistical model to calculate the 
# significance of the difference between the treatments and needs metadata 
# to define the treatments.

# The meta data needed for the frog data was information about which columns 
# were in which treatment group and which sibling group and we had that 
# information in a file. Similarly, here we need information on which columns 
# are from which cell type. Instead of having this is a file, we will create 
# a vector that indicates which column belongs to which cell type.

# Create a vector that indicates which column belongs to which cell type:

cell_type <- rep(c("prog","hspc"), 
                 times = c(length(prog) - 1,
                           length(hspc) - 1))
# The number of times each cell type is repeated is the number of columns 
# in that cell type minus 1. This is because we have removed the column with 
# the gene ids. Do check that the length of the `cell_type` vector is the 
# same as the number of columns in the `prog_hspc` dataframe.

# Run the differential expression analysis:
res_prog_hspc <- findMarkers(prog_hspc, 
                             cell_type)


# The dataframe `res_prog_hspc$prog` is log prog - log hspc (i.e.,Prog/HSPC). 
# This means 
# -   Positive fold change: prog is higher than hspc
# -   Negative fold change: hspc is higher than prog
# 
# The dataframe `res_prog_hspc$hspc` is log hspc - log prog (i.e., HSPC/Prog). 
# This means 
# -   Positive fold change: hspc is higher than prog
# -   Negative fold change: prog is higher than hspc


# Write the results to file:
data.frame(res_prog_hspc$prog, 
           ensembl_gene_id = row.names(res_prog_hspc$prog)) |> 
  write_csv("omics/week-4/results/prog_hspc_results.csv")






# OMICS 3 VISUALISE AND INTERPRET -----------------------------------------

library(tidyverse)
library(conflicted)
conflict_prefer("filter", "dplyr")


# import the normalised counts
prog <- read_csv("omics/week-5/data-raw/surfaceome_prog.csv")
hspc <- read_csv("omics/week-5/data-raw/surfaceome_hspc.csv")

# combine into one dataframe dropping one of the gene id columns
prog_hspc <- bind_cols(prog, hspc[-1])

# import the DE results
prog_hspc_results <- read_csv("omics/week-5/results/prog_hspc_results.csv")


# merge stats results with normalise values
prog_hspc_results <- prog_hspc_results |> 
  left_join(prog_hspc, by = "ensembl_gene_id")

# LINKING INFORMATION -----------------------------------------------------
# add the gene information using biomart

# [Ensembl](https://www.ensembl.org/index.html) [@birney2004] is a 
# bioinformatics project to organise all the biological information around
# the sequences of large genomes. The are a large number of databases 
# but [BioMart](https://www.ensembl.org/info/data/biomart/index.html) 
# [@smedley2009] provides a consistent interface to the material. 
# There are web-based tools to use these but the R package **`biomaRtr`**
# [@biomaRt] enables you to rapidly access and integrate information 
# into your R data structures.

library(biomaRt)

# Connect to the mouse database

ensembl <- useMart(biomart = "ensembl", 
                   dataset = "mmusculus_gene_ensembl")


# See what info we can retrieve:
listAttributes(mart = ensembl) |> View()

# We need
# ensembl_gene_id: because we will need to merge the info
# external_gene_name: name for gene
# description: description	

gene_info <- getBM(filters = "ensembl_gene_id",
                   attributes = c("ensembl_gene_id",
                                  "external_gene_name",
                                  "description"),
                   values = prog_hspc_results$ensembl_gene_id,
                   mart = ensembl)

# Notice the dataframe returned only has 279 rows, not 280. Which one is missing?

prog_hspc_results |> dplyr::select(ensembl_gene_id) |> 
  dplyr::filter(!ensembl_gene_id %in% gene_info$ensembl_gene_id)

# Error:
#   ! [conflicted] select found in 2 packages.
# Either pick the one you want with `::`:
# • biomaRt::select
# • dplyr::select
# Or declare a preference with `conflicts_prefer()`:
# • conflicts_prefer(biomaRt::select)
# • conflicts_prefer(dplyr::select)
# Run `rlang::last_trace()` to see where the error occurred.

prog_hspc_results |> dplyr::select(ensembl_gene_id) |> 
  filter(!ensembl_gene_id %in% gene_info$ensembl_gene_id)

# We might want to look that up. Google it. Let's worry about it 
# later if it turns out to be something important.

# merge the gene info with the results
prog_hspc_results <- prog_hspc_results |> 
  left_join(gene_info, by = "ensembl_gene_id")


# We now have datatframe with all the info we need, normalised counts, 
# log2 normalised counts, statistical comparisons with fold changes and p values, 
# information about the gene
# other than just the id


# save the most sig genes to file
prog_hspc_results_sig0.01 <- prog_hspc_results |> 
  filter(FDR <= 0.01)
# 168 genes

# write to csv file
write_csv(prog_hspc_results_sig0.01, 
          file = "omics/week-5/results/prog_hspc_results_sig0.01.csv")

prog_hspc_results_sig0.05 <- prog_hspc_results |> 
  filter(FDR <= 0.05)

# 182 genes

# write to csv file
write_csv(prog_hspc_results_sig0.05, 
          file = "omics/week-5/results/prog_hspc_results_sig0.05.csv")


# pca
# we do this on the log2 transformed normalised counts or the regularized the
# log transformed counts

# transpose the data
# we are reducing he number of dimensions from 280
prog_hspc_trans <- prog_hspc_results |> 
  dplyr::select(starts_with(c("Prog_", "HSPC_"))) |>
  t() |> 
  data.frame()

colnames(prog_hspc_trans) <- prog_hspc_results$ensembl_gene_id

# just for indep study before
# prog_hspc_trans$cell_id <- row.names(prog_hspc_trans)
# prog_hspc_trans <- prog_hspc_trans |> 
#   extract(cell_id, 
#           remove = FALSE,
#           c("cell_type", "cell_number"),
#           "([a-zA-Z]{4})_([0-9]{3})")
# 
# fig <- prog_hspc_trans |> ggplot(aes(x = ENSMUSG00000028639,
#                              y = ENSMUSG00000024053, colour = cell_type)) +
#   geom_point() +
#   # geom_text(aes(label = cell_id),
#   #           vjust = -1, size = 3) +
#   scale_x_continuous(expand = c(0.05,0.05)) +
#   scale_y_continuous(expand = c(0.05,0.05)) +
#   theme_classic() +
# theme(legend.position = "none")
# 
# 
# ggsave("omics/week-5/images/why_pca_mouse.png",
#        plot = fig,
#        width = 4, height = 4)


# perform PCA using standard functions
pca <- prog_hspc_trans |>
  prcomp(scale. = TRUE) 

summary(pca)
# Importance of components:
#                        PC1     PC2     PC3     PC4     PC5     PC6     PC7     PC8
# Standard deviation     4.3892 3.08797 2.25263 2.13943 1.96659 1.76697 1.62753 1.47668
# Proportion of Variance 0.0688 0.03406 0.01812 0.01635 0.01381 0.01115 0.00946 0.00779
# Cumulative Proportion  0.0688 0.10286 0.12098 0.13733 0.15114 0.16229 0.17175 0.17954

pca_labelled <- data.frame(pca$x,
                           cell_id = row.names(prog_hspc_trans))

# add the cell type information
# so we can label points 
# split cell_id into cell type and replicate and keep cell_id column

pca_labelled <- pca_labelled |> 
  extract(cell_id, 
          remove = FALSE,
          c("cell_type", "cell_number"),
          "([a-zA-Z]{4})_([0-9]{3})")



pca <- pca_labelled |> 
  ggplot(aes(x = PC1, y = PC2, 
             colour = cell_type)) +
  geom_point(alpha = 0.4) +
  scale_colour_viridis_d(end = 0.8, begin = 0.15,
                         name = "Cell type") +
  theme_classic()

# Fairly good separation of cell types but plenty of overlap

ggsave("omics/week-5/figures/prog_hspc-pca.png",
       plot = pca,
       height = 3, 
       width = 4,
       units = "in",
       device = "png")


# tSNE ??
# heatmap
library(heatmaply)
# we will use the most significant genes
# on a random subset of the cells
mat <- prog_hspc_results_sig0.01 |> 
  dplyr::select(starts_with(c("Prog", "HSPC"))) |>
  dplyr::select(sample(1:1499, size = 70)) |>
  as.matrix()

rownames(mat) <- prog_hspc_results_sig0.01$external_gene_name

n_cell_clusters <- 2
n_gene_clusters <- 2
heatmaply(mat, 
          scale = "row",
          hide_colorbar = TRUE,
          k_col = n_cell_clusters,
          k_row = n_gene_clusters,
          label_names = c("Gene", "Cell id", "Expression (normalised, log2)"),
          fontsize_row = 7, fontsize_col = 10,
          labCol = colnames(mat),
          labRow = rownames(mat),
          heatmap_layers = theme(axis.line = element_blank()))
# will take a few mins to run, and longer to appear in the viewer
# separation is not as strong as for the frog data


# volcano plot


# colour the points if FDR < 0.05
# and prog_hspc_results > 1
library(ggrepel)


prog_hspc_results <- prog_hspc_results |> 
  mutate(log10_FDR = -log10(FDR),
         sig = FDR < 0.05,
         bigfc = abs(summary.logFC) >= 2) 

vol <- prog_hspc_results |> 
  ggplot(aes(x = summary.logFC, 
             y = log10_FDR, 
             colour = interaction(sig, bigfc))) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), 
             linetype = "dashed") +
  geom_vline(xintercept = 1, 
             linetype = "dashed") +
  geom_vline(xintercept = -1, 
             linetype = "dashed") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_colour_manual(values = c("gray",
                                 "pink",
                                 "deeppink")) +
  geom_text_repel(data = subset(prog_hspc_results, 
                                bigfc & sig),
                  aes(label = external_gene_name),
                  size = 3,
                  max.overlaps = 50) +
  theme_classic() +
  theme(legend.position = "none")

ggsave("omics/week-5/figures/prog-hspc-volcano.png",
       plot = vol,
       height = 4.5, 
       width = 4.5,
       units = "in",
       device = "png")


# # just for the independent study slides
# vol <- prog_hspc_results |> 
#   ggplot(aes(x = summary.logFC, 
#              y = FDR)) +
#   geom_point() +
#   theme_classic() 
# ggsave("omics/week-5/images/volcano-why.png",
#        plot = vol,
#        height = 4.5, 
#        width = 4.5,
#        units = "in",
#        device = "png")