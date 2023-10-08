# OMICS 1 HELLO DATA ------------------------------------------------------

# Load tidyverse 
library(tidyverse)


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




# OMICS 3 VISUALISE AND INTERPRET -----------------------------------------