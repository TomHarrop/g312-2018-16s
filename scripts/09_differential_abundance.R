#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(phyloseq)

count_table_file <- "output/filter/count_table.txt"
taxonomy_file <- "output/taxonomy/aligned_otus.seed_v132.wang.taxonomy"

# make the count matrix
count_table <- fread(count_table_file)
count_matrix <- as.matrix(data.frame(count_table, row.names = "otu_id"))
otu <- otu_table(count_matrix, taxa_are_rows = TRUE)

# generate the taxonomy matrix
SplitTaxonomy <- function(x) {
    no_nums <- gsub("\\([[:digit:]]+\\)", "", x)
    strsplit(no_nums, ";")
}
head(fread(taxonomy_file, ))

taxonomy_table <- fread(taxonomy_file,
                        sep = '\t',
                        header = FALSE,
                        col.names = c("otu_id", "taxonomy"))
taxonomy_by_otu <- taxonomy_table[
    , transpose(SplitTaxonomy(taxonomy)),
    by = otu_id]
tax_levels <- c("Domain",
                "Phylum",
                "Class",
                "Order",
                "Family",
                "Genus")
setnames(taxonomy_by_otu,
         paste0("V", c(1:6)),
         tax_levels)
tax_matrix <- as.matrix(data.frame(taxonomy_by_otu, row.names = "otu_id"))
tax <- tax_table(tax_matrix)

# generate sample data
sample_table <- data.table(samplename = colnames(count_matrix))
sample_table[, c("population", "individual") := tstrsplit(samplename, "_")]
sd <- sample_data(data.frame(sample_table, row.names = "samplename"))
sd$samplename <- rownames(sd)

# generate phyloseq object
physeq <- phyloseq(otu, tax, sd)
physeq <- subset_taxa(physeq, !grepl("unclassified", Family)) 

# convert
dds <- phyloseq_to_deseq2(physeq, ~ population)

# run deseq
dds_wald <- DESeq(dds)
res_wald <- results(dds_wald,
                    contrast = c("population", "Lincoln", "Invermay"),
                    lfcThreshold = log(1.5, 2),
                    alpha = 0.1)

res_ordered <- res_wald[order(res_wald$padj), ]
tax_ordered <- tax[rownames(res_ordered)]

res_with_tax <- cbind(res_ordered, tax_ordered)

fwrite(data.table(data.frame(res_with_tax), keep.rownames = TRUE),
       "output/phyloseq/res_with_tax.csv")

