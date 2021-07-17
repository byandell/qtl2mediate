med_test <- function(med_ls, geno_max, phe_mx, kinship, cov_tar,
                     pos_Mbp, data_type, driver_med = NULL) {

  intermediate::mediation_test(target = phe_mx,
                               mediator = med_ls[[1]],
                               driver = geno_max,
                               annotation = med_ls[[2]],
                               covar_tar = cov_tar,
                               covar_med = med_ls$covar,
                               kinship = kinship,
                               driver_med = driver_med,
                               test = "wilc",
                               pos = pos_Mbp,
                               data_type = data_type)
}
misc_test_code <- function() {
  # see qtl2mediate.Rmd
  pos_Mbp <- mean(range(probs_obj$map))

  peak_mar <- qtl2::find_marker(probs_obj$map, chr_id, pos_Mbp)
  geno_max <- subset(probs_obj$probs, chr = chr_id, mar = peak_mar)[[1]][,,1]

  cov_tar <- qtl2shiny:::get_covar(
    covar, 
    dplyr::filter(analyses_tbl, pheno == pheno_name))
  driver_med <- probs_obj$probs[[chr_id]]

  phe_mx <- pheno_data[, pheno_name, drop = FALSE]
  
}

med_triad <- function(med_ls, geno_max, phe_mx, kinship, cov_tar, sdps, 
                     pattern, med_name, medID, haplos) {
  
  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  # Could use qtl2pattern::pull_mediator here.
  
  sdp <- sdps[qtl2pattern::sdp_to_pattern(sdps, haplos) == pattern]
  id <- med_ls[[2]]$id[med_ls[[2]][[medID]] == med_name]
  if(length(id) != 1)
    return(NULL)
  
  cov_tar <- covar_df_mx(cov_tar)
  cov_med <- covar_df_mx(med_ls$covar)
  
  intermediate::mediation_triad(target = phe_mx,
                              mediator = med_ls[[1]][, id, drop = FALSE],
                              driver = geno_max, 
                              covar_tar = cov_tar,
                              covar_med = cov_med,
                              kinship = kinship,
                              sdp = sdp, allele = TRUE)
}
misc_triad_code <- function() {
  # see qtl2mediate.Rmd
  snpinfo <- query_variants(chr_id, scan_window[1], scan_window[2])
  snpprobs <- qtl2shiny:::get_snpprobs(chr_id, peak_Mbp, window_Mbp,
                                       pheno_name, 
                                       probs_obj$probs,
                                       probs_obj$map,
                                       snpinfo)
  
  snp_scan_obj <- qtl2shiny:::scan1_covar(pheno_data, cov_tar,
                                          snpprobs$snpprobs, kinship,
                                          analyses_tbl,
                                          sex_type = "all")
  
  
  patterns <- summary(top_snps_tbl)
  haplos <- LETTERS[1:8]
  
  triad <- unique(med_test$best$triad)
  
  # Pick top trait with causal call.
  med_name <- dplyr::filter(med_test$best, .data$triad == "causal")[["id"]][1]
  pattern <- dplyr::filter(patterns, .data$pheno == med_name)$pattern[1]
  
  sdps <- unique(dplyr::filter(patterns, .data$pheno == med_name)$sdp)
  
  sdp <- sdps[qtl2pattern::sdp_to_pattern(sdps, haplos) == pattern]
  id <- med_ls[[2]]$id[med_ls[[2]][["id"]] == med_name]
  
  cov_tar <- qtl2shiny:::covar_df_mx(cov_tar)
  cov_med <- qtl2shiny:::covar_df_mx(med_ls$covar)
}
