#' @export
#' @importFrom qtl2 decomp_kinship fit1
#' @importFrom stringr str_replace
#' 
fitQtl2 <- function(driver,
                    target,
                    kinship = NULL,
                    ...) {

  if(is.null(rownames(driver)))
    rownames(driver) <- seq_len(nrow(driver))

  if(is.list(kinship))
    kinship <- kinship[[1]]
  
  
  out <- qtl2::fit1(driver, target, kinship, ...)
  
  # Replace lod names with LR
  names(out) <- stringr::str_replace(names(out), "lod", "LR")
  names(out) <- stringr::str_replace(names(out), "_LR", "LR")

  # Rescael to make them likelihoods (or likelihood ratios)
  out$LR <- out$LR * log(10)
  out$indLR <- out$indLR * log(10)
  
  # Add df for later use
  out$df <- ncol(driver) - 1
  
  # Residuals
  fitted <- rep(NA, length(target))
  names(fitted) <- if(is.matrix(target)) {
    rownames(target)
  } else {
    names(target)
  }
  fitted[names(out$fitted)] <- out$fitted
  out$resid <- target - fitted
  
  out
}