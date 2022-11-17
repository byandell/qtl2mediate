#' Locate triad for mediation
#' 
#' Given chr, pos and scan_window, target name and snpinfo,
#' find sdp for best pattern.
#' 
#' @param chr,pos,scan_window locators for region
#' @param pheno_name,pheno_data name and values of phenotype data
#' @param covar_tar optional covariates for target
#' @param genoprobs genotype probabilities
#' @param map map of genome
#' @param analyses_tbl analyses table
#' @param kinship kinship information (optional)
#' @param sex_type sex type: one of `c("A","I","F","M","all")`
#' @param snpinfo SNP information data frame (from user supplied \code{query_variants})
#' @param ... additional arguments
#' 
#' @export
#' @importFrom dplyr arrange desc
#' @importFrom rlang .data
#' 
triad_sdp <- function(chr, pos, scan_window,
                      pheno_name, pheno_data, covar_tar,
                      genoprobs, map, analyses_tbl, 
                      kinship = NULL,
                      sex_type = "I",
                      snpinfo = query_variants(chr, scan_window[1], scan_window[2]),
                      ...){
  
  if(!exists("query_variants")) {
    query_variants <- function(...) {NULL}
  }
  if(is.null(snpinfo)) {
    return(NULL)
  }
  # Compute SNP probabilities and update SNP info.
  snpprobs <- get_snpprobs(chr, pos, diff(scan_window) / 2,
                                         pheno_name, 
                                         genoprobs,
                                         map,
                                         snpinfo)
  
  # Scan SNPs based on 
  snp_scan_obj <- scan1covar(pheno_data[, pheno_name, drop = FALSE],
                                           covar_tar,
                                           snpprobs$snpprobs, kinship,
                                           dplyr::filter(analyses_tbl,
                                                         .data$pheno == pheno_name),
                                           sex_type = sex_type)
  
  # Can I use top_snps instead of top_snps_pattern?
  dplyr::arrange(
    qtl2::top_snps(snp_scan_obj, 
                   snpprobs$snpinfo,
                   lodcolumn = 1,
                   show_all_snps = FALSE),
    dplyr::desc(.data$lod))$sdp[1]
}

# Dummy routines. See qtl2::create_variant_query_func.
# query_variants <- function(...) {NULL}