#!/usr/bin/env R

library(tidyverse)
library(Maaslin2)
library(lubridate)
library(optparse)


option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="Input file path", metavar="character"),
  make_option(c("-m", "--metadata"), type="character", default=NULL,
              help="Metadata file path", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output dir path", metavar="character"),
  make_option(c("-f", "--fixed"), type="character", default=NULL,
              help="Fixed effects", metavar="character"),
  make_option(c("-d", "--dgroup"), type="character", default=NULL,
              help="Group variable", metavar="character"),
    make_option(c("-t", "--thread"), type="integer", default=1,
              help="Number of threads", metavar="integer")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)
cat("Reading metadata and input data...")
metapath <- opt$metadata
inputpath <- opt$input
metadata <- read.table(
    metapath,
    header = T,
    sep = "\t",
    row.names = 1
    )
abdf <- read.table(
    inputpath,
    header = T,
    sep = "\t",
    row.names = 1
  )
threads <- opt$thread
fixed_effects <- unlist(strsplit(opt$fixed, split = ","))
outdir <- opt$output
diet_group <- opt$dgroup
if (!(diet_group == "ALL")) {
  metadata <- metadata %>%
    dplyr::filter(Diet == diet_group)
  abdf <- abdf %>%
    rownames_to_column("sample") %>%                      # make sample column
    semi_join(metadata %>% rownames_to_column("sample"),  # keep only samples present in meta_sub
              by = "sample") %>%
    column_to_rownames("sample")   
}

diet_group_to_analyze <- unique(metadata[, "Diet"])
cat("Preparing to run Maaslin2 analysis... no sclaing\n")
cat(paste("Starting Maaslin2 analysis...on ", inputpath, "\n"))
cat(paste("Diet group: ", diet_group, "\n"))
cat(paste("Diet groups to analyze: ", paste(diet_group_to_analyze, collapse = ", "), "\n"))
cat(paste("number of samples : ", nrow(abdf), "\n"))
cat(paste("number of features : ", ncol(abdf), "\n"))
cat(paste("numberof matched_metadata samples : ", nrow(metadata), "\n"))
cat(paste("Fixed effects: ", paste(fixed_effects, collapse=", "), "\n"))


res <- Maaslin2(
  input_data =abdf,
  input_metadata = metadata,
  output = outdir,
  fixed_effects = fixed_effects,
  random_effects = c("mouse.ID", "Cohort", "HID", "ext.batch"),
  cores = threads,
  min_abundance = 0, min_prevalence = 0,
  normalization = "NONE", transform="NONE",
  plot_heatmap = FALSE,
  plot_scatter = FALSE)

cat(paste("Maaslin2 analysis completed. Results saved to", outdir, "\n"))



