# Load pckages ------------------------------------------------------------ 
library(tidyverse)
library(scran)
library(conflicted)

# OMICS 1 HELLO DATA ------------------------------------------------------

# HSPC --------------------------------------------------------------------


# Import ------------------------------------------------------------------

# 🐭 import the hspc data
hspc <- read_csv("data-raw/surfaceome_hspc.csv")

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
          file = "data-processed/hspc_summary_samp.csv")


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
          file = "data-processed/hspc_summary_gene.csv")



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
prog <- read_csv("data-raw/surfaceome_prog.csv")

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
          file = "data-processed/prog_summary_samp.csv")


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
          file = "data-processed/prog_summary_gene.csv")




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
  filter(sum(c_across(Prog_001:HSPC_852)) == 0)
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
# Find the genes that are 0 in every column for the Prog cells:
  
prog_hspc |> 
  rowwise() |> 
  filter(sum(c_across(HSPC_001:HSPC_852)) == 0)

# Note that if we knew there were some rows that were all zero across both 
# cell types, we would need to add 
# |> filter(sum(c_across(Prog_001:Prog_852)) != 0)

# Genes that are 0 in every column for the HSPC cells:
prog_hspc |> 
  rowwise() |> 
  filter(sum(c_across(HSPC_001:HSPC_852)) == 0)
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
  write_csv("results/prog_hspc_results.csv")






# OMICS 3 VISUALISE AND INTERPRET -----------------------------------------

