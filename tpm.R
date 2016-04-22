# a simple R function to calculate tpm (transcripts per million) for RNA-Seq data

TPM <- function (lengthkb, gene_count) {
  rpk <- gene_count/lengthkb
  total_rpk <- sum(rpk)
  tens_of_rpk<-total_rpk/1000000
  tpm<-rpk/tens_of_rpk
  return(tpm)
}
