#misc_test_code <- function() {
#  pos_Mbp <- mean(range(probs_obj$map))
#
#  peak_mar <- qtl2::find_marker(probs_obj$map, chr_id, pos_Mbp)
#  geno_max <- subset(probs_obj$probs, chr = chr_id, mar = peak_mar)[[1]][,,1]
#
#  cov_tar <- get_covar(
#    covar, 
#    dplyr::filter(analyses_tbl, pheno == pheno_name))
#  driver_med <- probs_obj$probs[[chr_id]]
#
#  phe_mx <- pheno_data[, pheno_name, drop = FALSE]
# 
#}

#misc_triad_code <- function() {
#  snpinfo <- query_variants(chr_id, scan_window[1], scan_window[2])
#  snpprobs <- get_snpprobs(
#    chr_id, peak_Mbp, window_Mbp,
#    pheno_name, 
#    probs_obj$probs,
#    probs_obj$map,
#    snpinfo)
#  
#  snp_scan_obj <- scan1_covar(
#    pheno_data, cov_tar,
#    snpprobs$snpprobs, kinship,
#    analyses_tbl,
#    sex_type = "all")
#  
#  patterns <- summary(top_snps_tbl)
#  haplos <- LETTERS[1:8]
#  
#  triad <- unique(med_test$best$triad)
#  
#  # Pick top trait with causal call.
#  med_name <- dplyr::filter(med_test$best, .data$triad == "causal")[["id"]][1]
#  pattern <- dplyr::filter(patterns, .data$pheno == med_name)$pattern[1]
#  
#  sdps <- unique(dplyr::filter(patterns, .data$pheno == med_name)$sdp)
#  
#  sdp <- sdps[sdp_to_pattern(sdps, haplos) == pattern]
#  id <- med_ls[[2]]$id[med_ls[[2]][["id"]] == med_name]
#  
#  cov_tar <- covar_matrix(cov_tar)
#  cov_med <- covar_matrix(med_ls$covar)
#}
