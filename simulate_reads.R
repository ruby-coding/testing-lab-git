#' @param n_snps The length of the sequence
#' @param n_true_seqs The number of sequences
#' @return A matrix of 0's and 1's with number of rows equal to n_snps and number of columns equal to n_true_seqs
make_true_seqs = function(n_snps, n_true_seqs, n_nonzero) {
  true_seqs = matrix(NA, nrow = n_true_seqs, ncol = n_snps)
  for (i in 1:n_snps) {
    snp_vals = rep(0, n_true_seqs)
    snp_vals[sample(n_true_seqs, n_nonzero)] = 1
    true_seqs[, i] = snp_vals
  }
  return(true_seqs)
}