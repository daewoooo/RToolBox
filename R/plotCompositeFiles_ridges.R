## This script will plot composite files in a circular layput ##
################################################################
## Load require libraries
library(primatR)
library(dplyr)
library(ggbio)
library(biovizBase)
library(BSgenome.Hsapiens.UCSC.hg38)
library(VariantAnnotation)

## Load SD annotation for SNV blacklisting
SDs.df <- read.table(file.path("/media/porubsky/DavStore/Diploid_assembly/GRCh38_annot/segDup_GRCh38.gz"))
SDs.df <- SDs.df[,c(2,3,4,27)]
colnames(SDs.df) <- c('seqnames', 'start', 'end', 'fracMatch')
SDs.gr <- makeGRangesFromDataFrame(SDs.df)
SDs.gr <- reduce(SDs.gr)
SDs.gr$ID <- 'SDs'

## Load HGSVC2 sample description
samples <- read.table("/media/porubsky/DavStore/HGSVC/HGSVC2_samples_populations.csv", sep = ',', header = FALSE, stringsAsFactors = FALSE)
colnames(samples) <- c('sample', 'pop', 'superpop', 'sex', 'source')
samples$id.num <- gsub("\\D", "", samples$sample)

## Load 1000G colors
sample.col <- read.table("/media/porubsky/DavStore/HGSVC/Plot_colors/U24_Subpopulations_Colors.csv", sep = ',', header=TRUE, stringsAsFactors = FALSE, comment.char = '&')
samples$pop.col <- sample.col$pop_color[match(samples$pop, sample.col$Population)]
samples$superpop.col <- sample.col$Spop_color[match(samples$superpop, sample.col$SuperPopulation)]

## Plot number of reads in each composite file
composite.files <- list.files(path = "/media/porubsky/DavStore/HGSVC/HGSVC2_U24_breakpoints_th.2_minR10_back.02", pattern = "_th.2_minR10_back.02.RData", full.names = TRUE)
#composite.files <- c("/media/porubsky/DavStore/Diploid_assembly/HG00733_toGRCh38/BreakpointR_results/syncReads_HG00733.RData", composite.files[1])
chromosomes <- paste0('chr', c(1:22, 'X'))

## Load non-redundant set of inversions
highlight.gr <- get(load("/media/porubsky/DavStore/HGSVC/INVcalls_final/nonred_inversions_n32.RData"))
highlight.gr  <- highlight.gr[,0]
highlight.gr$ID <- 'INV'

## Merge annotation ranges
annot.gr <- c(SDs.gr, highlight.gr)

## Plot number of reads in each composite file
chromosomes <- paste0('chr', c(1:22, 'X'))
region <- GRanges(seqnames='chr16', ranges=IRanges(start=12742805, end=37686164))
bsgenome <- BSgenome.Hsapiens.UCSC.hg38

## Get composite files to plot
composite.files <- list.files(path = "/media/porubsky/DavStore/HGSVC/HGSVC2_Compfiles/", pattern = "\\.RData", full.names = TRUE)

## Get sample names from VCFs
vcfFiles <- list.files("/media/porubsky/DavStore/HGSVC/PAV_SNVcallset/", pattern = '\\.vcf$', full.names = TRUE)
vcf.samples <- sapply(basename(vcfFiles), function(x) strsplit(x, '_')[[1]][1])
## Match IDs with VCF samples
composite.files.filt <- list()
for (i in seq_along(composite.files)) {
  composite.file <- composite.files[i]
  filename <- basename(composite.file)
  sample.id <- gsub(filename, pattern = "syncFrags_|\\.RData", replacement = '')
  sample.numeric.id <- gsub(sample.id, pattern = "[A-Z]", replacement = '', ignore.case = TRUE)
  sample.numeric.id <- gsub(sample.numeric.id, pattern = '002', replacement = '24385')
  if (any(grepl(vcf.samples, pattern = sample.numeric.id))) {
    composite.files.filt[[length(composite.files.filt) + 1]] <- composite.file
  }
}
composite.files.filt <- unlist(composite.files.filt, use.names = FALSE)

composite.files <- composite.files.filt  
reads.grl <- GRangesList()
for (i in seq_along(composite.files)) {
  composite.file <- composite.files[i]
  filename <- basename(composite.file)
  sample <- gsub(filename, pattern = 'syncFrags_|\\.RData', replacement = '')
  ## Filter by vcf.samples
  sample.numeric.id <- gsub(sample, pattern = "[A-Z]", replacement = '', ignore.case = TRUE)
  if (sample.numeric.id == '002') {
    sample.numeric.id <- '24385'
  }
  sample.vcf <- vcf.samples[grep(vcf.samples, pattern = sample.numeric.id)]
  ## Load composite file reads
  directional.reads <- get(load(composite.file))
  ## Get reads from ROI
  roi.reads <- subsetByOverlaps(directional.reads, region)
  roi.reads$sample <- sample.vcf
  reads.grl[[i]] <- roi.reads
}
reads.gr <- unlist(reads.grl, use.names = FALSE)
plt.df <- as.data.frame(reads.gr)
plt.df$midpoint <- plt.df$start + ( (plt.df$end - plt.df$start) %/% 2 )

## Add population ID and color
plt.df$pop.col <- samples$pop.col[match(plt.df$sample, samples$sample)]
plt.df$superpop.col <- samples$superpop.col[match(plt.df$sample, samples$sample)]
## Add population IDs
plt.df$pop <- samples$pop[match(plt.df$sample, samples$sample)]
plt.df$superpop <- samples$superpop[match(plt.df$sample, samples$sample)]

## Sort by population
superpop.ord <- c('AFR', 'SAS', 'EAS', 'EUR', 'AMR')
plt.df$superpop <- factor(plt.df$superpop, levels = superpop.ord)

# plt <- ggplot(plt.df, aes(x = start, y = sample, fill = strand, height=..density..)) +
#     geom_density_ridges(stat = "density", alpha=0.75) +
#     #facet_grid(. ~ roi, scales = 'free') +
#     facet_grid(superpop ~ roi, scales = 'free') +
#     scale_x_continuous(expand = c(0,0)) +
#     scale_fill_manual(values = c("paleturquoise4", "sandybrown")) +
#     xlab("Genomic position (bp)") +
#     #theme(axis.text.x = element_text(colour = a)) +
#     theme_bw()

## Plot chromosome16 example
library(ggridges)
plt <- ggplot(plt.df, aes(x = start, y = sample, fill = strand)) +
  geom_density_ridges(stat = "binline", binwidth=50000, scale=6, alpha=0.75) +
  #facet_grid(sample ~ ., scales = 'free') +
  scale_x_continuous(expand = c(0,0), labels = comma) +
  scale_fill_manual(values = c("paleturquoise4", "sandybrown")) +
  xlab("Genomic position (bp)") +
  theme_bw()
## Export plot
destination <- "/media/porubsky/DavStore/HGSVC/CompositeFilesPlots/compositeFile_readDensity_chr16_ridges.png"
png(filename = destination, width = 10, height = 5, units = 'in', res = 300)
plt
dev.off()
