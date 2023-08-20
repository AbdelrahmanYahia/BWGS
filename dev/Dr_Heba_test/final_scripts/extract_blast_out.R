defaultW <- getOption("warn") 
options(warn = -1) 

suppressMessages(library("tidyverse", quietly = T))
suppressMessages(library("optparse", quietly = T))

#######################################
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="blast out file of frmt 6", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL, 
              help="output table (.csv)", metavar="character"),
  make_option(c("--identity"), type="numeric", default=97, 
              help="minimum identity", metavar="numeric"),
  make_option(c("--coverage"), type="numeric", default=90, 
              help="minimum coverage", metavar="numeric"),
  make_option(c("--evalue"), type="numeric", default=1e-4, 
              help="minimum e-value", metavar="numeric")
  );

#######################################
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
write(paste0("\033[0;", 32, "m","Generating Filterred table...","\033[0m"), stderr())
#######################################
split_col <- function(pattern, col, part){ # function for splitting a string into a vector according to a certain pattern
  col.chr <- as.character(col)
  lst <- sapply(strsplit(col.chr, split = pattern),`[`, part)
  vctr <- unlist(lst)
  return(vctr)
}

data <- read.table(opt$input, header=F)
# add column names
outfmt <- "qseqid sseqid pident qlen slen sseq length mismatch gaps qstart qend sstart send evalue bitscore"
col_names <- unlist(strsplit(outfmt,split = " "))
colnames(data) <- col_names

plas_out_flt <- data%>% # select rows wih perecnt-identity more than or equal to 97
  filter(pident >= opt$identity & evalue <= opt$evalue)


# innitiate empty vectors 
{
  ttl_smpl_rds <- nrow(plas_out_flt)
  unique_plas<- unique(plas_out_flt[,c("sseqid","slen")])
  nrow_unique <- nrow(unique_plas)
  plasmids <- numeric(nrow_unique)
  no_reads <- numeric(nrow_unique)
  sample_id <- numeric(nrow_unique)
  coverage <- numeric(nrow_unique)
}

for (i in 1:nrow_unique){
  plasmid_id <- as.character(unique_plas$sseqid[i])
  rows<-which((plas_out_flt$sseqid)==plasmid_id) # indeces of the origin of replication
  plasmids[i] <- plasmid_id # add to vector
  positions <- vector() # empty vector
  seq_len <- unique_plas$slen[i] #  length of the sequence
  for (n in rows){
    start <- plas_out_flt$sstart[n] # start position of the read
    end <- plas_out_flt$send[n] # end position of the read
    read <- start:end # generate sequence of position between the sgtart and the end
    positions <- append(positions,read) # add the sequence of positions of aligned reads 
  }
  # add number of aligned reads to a vector of the number of aligned reads to each sequence
  no_reads[i] <- length(rows)
  # add number of aligned reads to a vector of the number of aligned reads to each sequence
  total_positions <- length(unique(positions))
  # calcualte the coverage of each sequence, then add it to a vector
  coverage[i] <- (total_positions/seq_len)*100
}

plasmid_id <- split_col("\\:", plasmids, 1) # extract the id of each plasmid
# combine all the vectors into a data frame
plas_ori_cov <- data.frame(plasmid_id,coverage,no_reads)

#filter out sequences with coverage less than 90
plas_result <- plas_ori_cov%>%
  filter(coverage > opt$coverage)

write.csv(plas_result, opt$output, quote = F)
