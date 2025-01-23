source('methylation_R_utils.R')

library(GenomicRanges)
library(readr)
library(bsseq)
library(stats)
library(tibble)
library(dplyr)

cyto = read.csv("cytobands_GRCh38", header=FALSE, sep = '\t')
cyto_window <-  GRanges(seqnames = cyto$V1, ranges = IRanges(start = cyto$V2, end = cyto$V3))
seqnames_char <- as.character(seqnames(cyto_window))


for (j in c((1:22),'X','Y')){

  chrom <- paste0('chr',j)
  
  print(chrom)
  
  c <- read.csv(paste0('loci.',chrom,'.bed'), header = FALSE, sep = '\t')
  window <- GRanges(seqnames = c$V1, ranges = IRanges(start = c$V2, end = c$V2))
  chrom_indices <- grepl(chrom, seqnames_char)
  cyto_chrom <- cyto_window[chrom_indices]
  
  for (i in 1:length(cyto_chrom)){
    d <- paste0(chrom,'/',chrom,'.cyto.',i,'.mfreq.txt')
    ov <- GenomicRanges::intersect(cyto_chrom[i], window)
    df_matches <- c[c$V2 %in% start(ov), ]
    
    if (isEmpty(df_matches)) next
    
    write.table(df_matches, row.names = FALSE, col.names = FALSE, sep = '\t',
                file = d, quote = FALSE)
    bsobj <- read.bismark(d)
    gpc <- BSmooth(bsobj, ns=10, h=100)
    gpc.cov <- getCoverage(gpc,type="Cov",what="perBase")[,1]
    
    keepi <- which(gpc.cov>5)
    gpc.keep <- gpc[keepi,]
    
    if (length(granges(gpc.keep))==0) next
    
    #regions <- gpcPeakCaller(gpc.keep)
    
    minwin = 50
    a = 0.01
    cutoff = NULL
    qcutoff = 0.99
    bsobj <- gpc.keep
    gpc.meth <- getMeth(bsobj, type="smooth", what="perBase")[,1]
    # remove NA
    keepi <- !is.na(gpc.meth)
    bsobj <- bsobj[keepi,]
    gpc.meth <- gpc.meth[keepi]
    gpc.cov <- getCoverage(bsobj,type="Cov",what="perBase")[,1]
    gpc.m <- getCoverage(bsobj,type="M",what="perBase")[,1]
    ## compare to baseline ----
    baseline <- median(gpc.meth,na.rm = T)
    gpc.diff <- gpc.meth - baseline
    if ( is.null(cutoff)) {
      cutoff <- quantile(gpc.diff,qcutoff,na.rm = T)
    }
    gpc.direction <- ifelse(gpc.diff > cutoff, 1, -1) # cut off by qcutoff
    gpc.gr <- granges(bsobj)
    chrs <- as.character(seqnames(gpc.gr))
    pos <- start(gpc.gr)
    ## find regions of + ----
    regions <- bsseq:::regionFinder3(gpc.direction,chrs,pos)$up
    regions <- as_tibble(regions)
    ## then add average and peak accesibility, along with coverage, then binomial test ---
    
    if (length(regions) == 0) next
    
    regions <- regions %>%
      rowwise() %>%
      mutate(
        coverage = sum(gpc.cov[idxStart:idxEnd]),
        methylated = sum(gpc.m[idxStart:idxEnd]),
        average = mean(gpc.meth[idxStart:idxEnd]),
        peak = max(gpc.meth[idxStart:idxEnd]),
        p.value = binom.test(methylated,coverage,baseline,alternative = "greater")$p.value
      ) %>%
      ungroup()
    ## multiple test adjusting using FDR
    regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")
    
    ## significance based on :
    ## 1. width
    ## 2. alpha
    ## 3. minimum peak height
    minpeak <- 0 #baseline + 2 * cutoff 
    regions <- regions %>%
      mutate(
        width = end - start + 1)
    regions <- regions %>%
      mutate(
        sig = ifelse(
          adjusted.pval <= a &
            width >= minwin &
            peak >= minpeak
          ,"sig","insig"))
    
    if (all(regions$peak <= 0.699999999999999)) next
    
    
    #regions %>%
    #  filter(sig == "sig", peak >= 0.7) %>%
    #  summarize( min(average), min(peak), n())
    regions %>%
      dplyr::filter(sig == "sig") 
    regions %>%
      dplyr::filter(peak >= 0.7) %>%
      summarize(min(average), min(peak), n())
    table(regions$sig)
    regions.sig <- regions %>%
      dplyr::filter(sig == "sig")
    table(regions$sig)
    regions.sig <- regions %>%
      dplyr::filter(sig == "sig")
    out <- regions.sig %>%
      dplyr::select(-cluster,-sig) %>%
      dplyr::rename(Chromosome = "chr",
                    Start = "start",
                    End = "end",
                    Num_Sites = "n",
                    Observations = "coverage",
                    Methylated = "methylated",
                    Average_Accessibility = "average",
                    Maximum_Accessibility = "peak",
                    Region_Width = "width"
      )
    
    
    outpath <- file.path(paste0(chrom,'/',chrom,'.',i,'.peak.tsv'))
    write_tsv(out,outpath)
    rds_path <- file.path(paste0(chrom,'/',chrom,'.',i,'.peaks.Rds'))
    saveRDS(regions, rds_path)
    
    data <- read.table(paste0(chrom,'/',chrom,'.',i,'.peak.tsv'), header = TRUE, sep = '\t')
    write.table(data, paste0(chrom,'.all.peaks.tsv'), sep="\t", row.names=FALSE, col.names=FALSE, append=TRUE, quote = FALSE)
    
    
  }
  
  
}

for (k in c((1:22),'X','Y')) {
  
  chrom <- paste0('chr',k)
  print(chrom)
  df <- read.table(paste0(chrom,'.all.peaks.tsv'), header = TRUE, sep = '\t')
  write.table(df, paste0(chrom,'.all.peaks.tsv'), sep="\t", row.names=FALSE, col.names=FALSE, quote = FALSE)
  
}

for (m in c((1:22),'X','Y')){
  chro <- paste0('chr',m)
  all <- read.table(paste0(chro,'.all.peaks.tsv'), header = TRUE, sep = '\t')
  write.table(all, paste0('all.peaks.tsv'), sep="\t", row.names=FALSE, col.names=FALSE,append = TRUE, quote = FALSE )
}

for (bb in c((1:22),'X','Y')){
  print(bb)
}
