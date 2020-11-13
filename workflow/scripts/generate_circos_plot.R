# R script to generate Circos plot based on BAM and fa.gz files

# Author: Maurits Evers


suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(circlize))


# Snakemake parameters
bam_file <- snakemake@input[["bam"]]
rDNA_file <- snakemake@input[["rDNA"]]
cyto_file <- snakemake@input[["cyto"]]
pdf_file <- snakemake@output[["pdf"]]
png_file <- snakemake@output[["png"]]


# Read data
bam <- BamFile(bam_file)
rDNA <- readDNAStringSet(rDNA_file)
cyto <- read_delim(
	cyto_file,
	col_types = "ciicc",
	delim = "\t",
	col_names = c("chrom", "start", "end", "name", "giemsa_stain"))

# Get chromosome/rDNA sequence lengths from BAM and rDNA sequencing files
df_size <- bam %>%
	scanBamHeader() %>%
	pluck("targets") %>%
	c(rDNA = width(rDNA)) %>%
	enframe()


# Generate cytoband `data.frame` based on UCSC data and add rDNA entry
df_cyto <- cyto %>%
	mutate(chrom = str_remove(chrom, "chr")) %>%
	mutate(chrom = if_else(chrom == "M", "MT", chrom)) %>%
	filter(chrom %in% df_size$name) %>%
	add_row(
		chrom = "rDNA", start = 1, end = width(rDNA),
		name = "", giemsa_stain = "") %>%
	as.data.frame()


# Generate rDNA to genome read map based on BAM reads
df_map <- bam %>%
	scanBam() %>%
	pluck(1) %>%
	keep(names(.) %in% c("qname", "rname", "pos", "qwidth")) %>%
	as_tibble() %>%
	separate(qname, c("tmp", "read_no", "from_pos"), sep = "_") %>%
	select(-tmp, -read_no) %>%
	mutate(from_chr = "rDNA", from_pos = str_remove(from_pos, "pos")) %>%
	separate(
		from_pos, into = c("from_start", "from_end"),
		sep = "-", convert = TRUE) %>%
	rename(to_chr = rname, to_start = pos, to_end = qwidth) %>%
	mutate(to_end = to_end + to_start) %>%
	select(from_chr, from_start, from_end, to_chr, to_start, to_end) %>%
	as.data.frame()


# Make sure that chromosome names from `df_size` match those from `df_cyto` and
# `df_map`
valid_chr <- intersect(df_size$name, df_cyto$chrom)
df_size <- df_size %>% filter(name %in% valid_chr)
df_cyto <- df_cyto %>% filter(chrom %in% valid_chr)
df_map <- df_map %>% filter(to_chr %in% valid_chr)


# Generate title from BAM input file
chunks <- bam_file %>%
	str_split("[/_]") %>%
	unlist() %>%
	str_remove_all("(len|step|\\.bam)")
title <- sprintf("Reference: %s (%s); rDNA fragments from sliding windows (length = %s bp, step size = %s bp)",
	str_to_title(chunks[3]), chunks[4], chunks[7], chunks[8])


# Save as PDF
pdf(pdf_file)
circos.initializeWithIdeogram(
	df_cyto,
	chromosome.index = valid_chr,
	sector.width = c(rep(1, length(valid_chr) - 1), 10))
circos.genomicLink(
	df_map[1:3], df_map[4:6],
	col = rgb(0, 0, 0, max = 255, alpha = 20), border = "black", lwd = 0.01)
title(title, cex.main = 0.65)
dev.off()


# Save as PNG
png(png_file, height = 1800, width = 1800, res = 300)
circos.initializeWithIdeogram(
	df_cyto,
	chromosome.index = valid_chr,
	sector.width = c(rep(1, length(valid_chr) - 1), 10))
circos.genomicLink(
	df_map[1:3], df_map[4:6],
	col = rgb(0, 0, 0, max = 255, alpha = 20), border = "black", lwd = 0.01)
title(title, cex.main = 0.65)
dev.off()
