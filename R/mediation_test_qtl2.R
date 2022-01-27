# Mediation tests for qtl2 data across index
#
#' Test mediation across set of indexed drivers for `qtl2` data
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param genoprobs genoprob object of class [qtl2::calc_genoprob()]
#' @param map list of map positions
#' @param chr,pos chromosome and position of driver
#' @param ... additional parameters for [mediation_test()]
#'
#' @importFrom qtl2 find_marker
#' @importFrom intermediate mediation_test
#' 
#' @return Object of class `mediation_qtl2`, which inherits from class [mediation_index()]
#' 
#' @examples 
#' \donttest{
#' dirpath <- "https://raw.githubusercontent.com/rqtl/qtl2data/master/DOex"
#' 
#' # Read DOex example cross from 'qtl2data'
#' DOex <- subset(qtl2::read_cross2(file.path(dirpath, "DOex.zip")), chr = "2")
#' 
#' # Download genotype probabilities
#' tmpfile <- tempfile()
#' download.file(file.path(dirpath, "DOex_genoprobs_2.rds"), tmpfile, quiet=TRUE)
#' pr <- readRDS(tmpfile)
#' unlink(tmpfile)
#' 
#' med_test <- qtl2mediate::mediation_test_qtl2(
#'   target = DOex$pheno,
#'   mediator = DOex$pheno,
#'   annotation = med_ls[[2]],
#'   covar_tar = cov_tar,
#'   covar_med = med_ls$covar,
#'   genoprobs = probs_obj$probs,
#'   map = probs_obj$map,
#'   chr = chr_id,
#'   pos = pos_Mbp,
#'   kinship = kinship)
#' }
#' @export
#'
mediation_test_qtl2 <- function(target,
                                mediator,
                                annotation,
                                covar_tar = NULL,
                                covar_med = NULL,
                                genoprobs,
                                map,
                                chr,
                                pos,
                                ...) {
  
  # Get driver as genoprobs based on chr and pos.
  peak_mar <- qtl2::find_marker(map, chr, pos)
  driver <- subset(genoprobs, chr = chr, mar = peak_mar)[[1]][,,1]
  
  # Get driver_med for use with annotation data on mediators.
  driver_med <- genoprobs[[chr]]
  
  intermediate::mediation_test(target,
                               mediator,
                               driver,
                               annotation,
                               covar_tar,
                               covar_med,
                               driver_med = driver_med,
                               fitFunction = fitQtl2,
                               ...)
}