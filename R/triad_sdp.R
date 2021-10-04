#' Locate triad for mediation
#' 
#' Given chr, pos and scan_window, target name and snpinfo,
#' find sdp for best pattern.
#' 
#' @param chr,pos,scan_window locators for region
#' @param pheno_name,pheno_data name and values of phenotype data
#' @param genoprobs genotype probabilities
#' @param map map of genome
#' @param analyses_tbl analyses table
#' @param kinship kinship information (optional)
#' @param snpinf SNP information data frame (from user supplied \code{query_variants})
#' @param ... additional arguments
#' 
#' @export
#' @importFrom dplyr arrange desc
#' 
triad_sdp <- function(chr, pos, scan_window,
                      pheno_name, pheno_data,
                      genoprobs, map, analyses_tbl, 
                      kinship = NULL,
                      snpinfo = query_variants(chr, scan_window[1], scan_window[2]),
                      ...){
  
  # Compute SNP probabilities and update SNP info.
  snpprobs <- qtl2mediate:::get_snpprobs(chr, pos, diff(scan_window) / 2,
                                         pheno_name, 
                                         genoprobs,
                                         map,
                                         snpinfo)
  
  # Scan SNPs based on 
  snp_scan_obj <- qtl2mediate:::scan1covar(pheno_data[, pheno_name, drop = FALSE],
                                           cov_tar,
                                           snpprobs$snpprobs, kinship,
                                           dplyr::filter(analyses_tbl,
                                                         pheno == pheno_name),
                                           sex_type = "all")
  
  # Can I use top_snps instead of top_snps_pattern?
  dplyr::arrange(
    qtl2::top_snps(snp_scan_obj, 
                   snpprobs$snpinfo,
                   lodcolumn = pheno_name,
                   show_all_snps = FALSE),
    dplyr::desc(.data$lod))$sdp[1]
}
