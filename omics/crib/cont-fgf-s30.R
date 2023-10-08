
# OMICS 1 HELLO DATA ------------------------------------------------------

# Load tidyverse 
library(tidyverse)


# Import ------------------------------------------------------------------

# 🐸 import the s30 data
s30 <- read_csv("data-raw/xlaevis_counts_S30.csv")


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




# OMICS 3 VISUALISE AND INTERPRET -----------------------------------------


