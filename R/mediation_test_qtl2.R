# Mediation tests for qtl2 data across index
#
#' Test mediation across set of indexed drivers for `qtl2` data
#'
#' @param target vector or 1-column matrix with target values
#' @param mediator vector or 1-column matrix with mediator values
#' @param annotation optional annotation data frame for mediators
#' @param covar_tar optional covariates for target
#' @param covar_med optional covariates for mediator
#' @param kinship optional kinship matrix among individuals
#' @param genoprobs genoprob object of class [qtl2::calc_genoprob()]
#' @param map list of map positions
#' @param drop_lod drop in `LOD` (= `LR/log(10)`) to define index set
#' @param min_lod minimum lod to consider (default \code{3})
#' @param query_variant function to query variant database
#' @param cores number of cores to use
#' @param target_scan optional object from [qtl2::scan1snps()] for target (created if missing)
#' @param ... additional parameters for [mediation_index()]
#'
#' @importFrom qtl2 find_marker
#' @importFrom intermediate mediation_test
#' 
#' @return Object of class `mediation_qtl2`, which inherits from class [mediation_index()]
#' 
#' @examples 
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
                                kinship = NULL,
                                ...) {
  
  peak_mar <- qtl2::find_marker(map, chr, pos)
  driver <- subset(genoprobs, chr = chr, mar = marker)[[1]][,,1]
  driver_med <- genoprobs[[chr]]
  
  intermediate::mediation_test(target,
                               mediator,
                               driver,
                               annotation,
                               covar_tar,
                               covar_med,
                               driver_med = driver_med,
                               kinship = kinship,
                               ...)
}