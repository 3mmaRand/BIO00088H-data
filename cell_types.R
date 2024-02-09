# 13.1 Cell types data set
# 13.1.1 ST1
cell_types <- as.matrix(read.delim("omics/all_cell_types.txt", sep =
                                     "\t", row.names = 1))


# 13.1.2 ST2
cell_types[1:5,1:5]
dim(cell_types)


# 13.1.3 ST3
colnames(cell_types)


# 13.1.4 ST4
table(cell_types)
cell_types[,"LTHSC"]


# 13.1.5 ST5
sort(colSums(cell_types)) |> sum()


# 13.1.6 ST6




# 13.1.7 ST7
cell_labels <- sub("_.*", "", rownames(cell_types))
table(cell_labels)
# HSPC LT.HSC   Prog 
# 852    216    852 
# 701    155    798 orginal


# 13.1.8 ST8
colSums(cell_types[grep("LT.HSC", rownames(cell_types)),])
colSums(cell_types[grep("HSPC", rownames(cell_types)),])
colSums(cell_types[grep("Prog", rownames(cell_types)),])


# 13.1.9 ST10
length(colnames(lthsc))
length(rownames(cell_types))
length(intersect(colnames(lthsc),rownames(cell_types)))
# 13.2 Combining the gene expression data with the cell type information
# Back to questions


cell_is_cmp <- clean_cell_types[,"CMP"]==1
cmp_cell_names <- rownames(clean_cell_types[cell_is_cmp,])
cmp_expression <- gene_counts[,cmp_cell_names]



# my code ---------------------------------------------------------
library(tidyverse)
# import each of the data sets
# 🐭 import the data
hspc <- read_csv("omics/week-5/data-raw/surfaceome_hspc.csv")
prog <- read_csv("omics/week-5/data-raw/surfaceome_prog.csv")

# what are the name of the cells we have data for
hspc_cells <- colnames(hspc[-1])
prog_cells <- colnames(prog[-1])
all_cells <- c(hspc_cells, prog_cells)

# import the extra cell typ info
cell_types <- read_table("omics/all_cell_types.txt")


# subset the extra cell info using only cells for for we have expression data
cell_types <- cell_types |> 
  filter(cell %in% all_cells)


# write that to file
write_csv(cell_types, "omics/er_cell_types.csv")

# 
