#/usr/bin/Rscript
library(Biostrings)
f_dir <- "/mnt/ceph/454_data/TideHunter_tests/jirka/Cameor_v2/analysis_petr/tmp_test"
fasta_files <- dir(f_dir, pattern = ".fasta")

for (i in fasta_files){
  s <- readDNAStringSet(paste(f_dir, i, sep = "/"))
  if (length(s) < 5){
    next
  }
  km <-  oligonucleotideFrequency(s, width = 7, as.prob = TRUE) + oligonucleotideFrequency(reverseComplement(s), width = 7, as.prob = TRUE)
  d <- dist(km)
  mds <- cmdscale(d)
  plot(mds, main = i)
  Sys.sleep(1)
  hc <- hclust(d)
}



