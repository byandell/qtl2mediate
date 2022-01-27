#' Create list with mediators in region
#' 
#' @param chr_id,scan_window chromosome and start and end value
#' @param covar covariate data frame
#' @param map list or vector of map positions 
#' @param peaks table of peaks
#' @param analyses table of analyses
#' @param pheno_data matrix of phenotype data
#' @param drivers number of drivers (1 or 2; default is 2)
#' 
#' @importFrom rlang .data
#' @importFrom dplyr filter
#' 
#' @export
comediator_region <- function(pheno_name, chr_id, scan_window, 
                              covar, analyses_tbl, peaks, 
                              qtls = 2, pmap, pheno_data) {
  
  # This is specific to CCmouse.
  peaks <- dplyr::filter(peaks,
                         .data$pheno != pheno_name,
                         !(.data$pheno_group == "Islet.mRNA"))
  
  # Filter peaks and analyses to region and drop pheno_name
  peaks_local <- dplyr::filter(peaks,
                               .data$chr == chr_id,
                               .data$pos >= scan_window[1],
                               .data$pos <= scan_window[2])
  analyses_df <- dplyr::filter(analyses_tbl,
                               .data$pheno %in% peaks_local$pheno,
                               .data$pheno != pheno_name)
  
  # Read the phenos we need.
  phenos <- analyses_df$pheno
  pheno_data <- pheno_data[, phenos, drop = FALSE]

  # Create comediator object.
  out <- pheno_region(chr_id, scan_window, covar, pmap,
    peaks, analyses_tbl, pheno_data, drivers = qtls)
  
  out
}
#' Get mediators of same type as pheno_name
#' 
#' @export
comediator_type <- function(comed, peaks, pheno_name, doThis) {
  if(doThis & !is.null(comed)) {
    phe_type <- (dplyr::filter(peaks, .data$pheno == pheno_name))$pheno_type[1]
    not_type <- comed$annot$id[comed$annot$biotype != phe_type]
    comed$comediators <- comed$comediators[, not_type]
  }

  comed
}
