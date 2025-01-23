library(GenomicRanges)
library(readr)
library(bsseq)
library(stats)
library(tibble)
library(dplyr)
library(EnsDb.Hsapiens.v75)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)
library(data.table)

get_strand_sensitive_motif <- function(chromosome, position, upstream, downstream, genome, reverse = FALSE) {
  start_pos <- position - upstream
  end_pos <- position + downstream
  sequence <- getSeq(genome, chromosome, start=start_pos, end=end_pos)
  if (reverse) {
    sequence <- reverseComplement(sequence)
  }
  return(sequence)
}

genome <- BSgenome.Hsapiens.UCSC.hg38
upstream <- 3
downstream <- 3

# Initialize an empty data.table
newdf <- data.table(
  chr = character(),
  pos = numeric(),
  strand = character(),
  meth = integer(),
  unmeth = integer(),
  context = character(),
  motif = character()
)

# Read the total number of lines in the file
total_lines <- as.numeric(system("wc -l < cov2cyt.NOMe.GpC_report.txt", intern = TRUE))
chunk_size <- 100000  # Adjust the chunk size based on your memory capacity

# Read and process the file in chunks
tab <- fread('cov2cyt.NOMe.GpC_report.txt', sep = '\t', header = FALSE, nrows = chunk_size, select = 1:5)
chunk_counter <- 1

while (nrow(tab) > 0) {
  # Split the data by strand
  tab_plus <- tab[tab$V3 == '+']
  tab_minus <- tab[tab$V3 == '-']
  
  # Process the '+' strand
  if (nrow(tab_plus) > 0) {
    motifs_plus <- mapply(get_strand_sensitive_motif, 
                          chromosome = tab_plus[, 1], 
                          position = tab_plus[, 2], 
                          MoreArgs = list(upstream = upstream, downstream = downstream, genome = genome, reverse = FALSE))
    
    true_motifs_plus <- sapply(motifs_plus, function(motif) as.character(motif[3:5]))
    
    valid_indices_plus <- which(true_motifs_plus %in% c('GCC', 'GCT', 'GCA'))
    if (length(valid_indices_plus) > 0) {
      new_rows_plus <- data.table(
        chr = tab_plus[valid_indices_plus, 1],
        pos = tab_plus[valid_indices_plus, 2],
        strand = '+',
        meth = tab_plus[valid_indices_plus, 4],
        unmeth = tab_plus[valid_indices_plus, 5],
        context = 'GC',
        motif = true_motifs_plus[valid_indices_plus]
      )
      setnames(new_rows_plus, old = names(new_rows_plus), new = c("chr", "pos", "strand", "meth", "unmeth", "context", "motif"))
      newdf <- rbind(newdf, new_rows_plus)
    }
  }
  
  # Process the '-' strand
  if (nrow(tab_minus) > 0) {
    motifs_minus <- mapply(get_strand_sensitive_motif, 
                           chromosome = tab_minus[, 1], 
                           position = tab_minus[, 2], 
                           MoreArgs = list(upstream = upstream, downstream = downstream, genome = genome, reverse = TRUE))
    
    true_motifs_minus <- sapply(motifs_minus, function(motif) as.character(motif[3:5]))
    
    valid_indices_minus <- which(true_motifs_minus %in% c('GCC', 'GCT', 'GCA'))
    if (length(valid_indices_minus) > 0) {
      new_rows_minus <- data.table(
        chr = tab_minus[valid_indices_minus, 1],
        pos = tab_minus[valid_indices_minus, 2],
        strand = '-',
        meth = tab_minus[valid_indices_minus, 4],
        unmeth = tab_minus[valid_indices_minus, 5],
        context = 'GC',
        motif = true_motifs_minus[valid_indices_minus]
      )
      setnames(new_rows_minus, old = names(new_rows_minus), new = c("chr", "pos", "strand", "meth", "unmeth", "context", "motif"))
      newdf <- rbind(newdf, new_rows_minus)
    }
  }
  
  
  }
  
  
  tab <- fread('cov2cyt.NOMe.GpC_report.txt', sep = '\t', header = FALSE, nrows = chunk_size, skip = nrow(newdf), select = 1:5)
}

# Sort the final data frame
setorder(newdf, chr, pos)

write.table(newdf, file = 'loci.bed', sep = '\t', quote = FALSE, row.names = FALSE, col.names = FALSE)
                                
# Test the function with specific inputs
args <- list(
  chromosome = "chr1",
  position = 10482,
  upstream = upstream,
  downstream = downstream,
  genome = genome,
  reverse = FALSE
)
motif <- get_strand_sensitive_motif(args$chromosome, args$position, args$upstream, args$downstream, args$genome, args$reverse)
print(motif)

args$reverse <- TRUE
motif <- get_strand_sensitive_motif(args$chromosome, args$position, args$upstream, args$downstream, args$genome, args$reverse)
print(motif)
