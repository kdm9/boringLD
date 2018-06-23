require(doFuture)
require(foreach)
require(dplyr)

rho2halfmax = function(r, n, maxd=100000) {
  # This is a horrible hack as I can't seem to solve the below for d, and do
  # f(r, rsq, n)(r, rsq(d=1)/2, n)  to find distance of halfmaximal analytically
  #d = as.numeric(1:maxd)
  d = ceiling(c(1:99, 10^((200:900)/100)))
  r2 = ((10+r*d)/((2+r*d)*(11+r*d))) *(1+((3+r*d)*(12+12*r*d+(r*d)^2))/(n*(2+r*d)*(11+r*d)))
  d[min(which(r2 < max(r2)/2))]
}

#' window_halfmax
#'
#' @param bcf Indexed BCF or VCF file
#' @param region Region of genome
#' @param minMAF Minimum SNP minor allle freq
#' @param maxMissing Maximum SNP missing data rate
#'
#' @return data.frame with one row and columns rho, halfmax, region, nsnp
#' @export
window_halfmax = function(bcf, region, minMAF, maxMissing, geno=NULL) {
  if (is.null(geno)) {
    geno = bcf_getGTandAD(bcf, region, minMAF = minMAF, maxMissing = maxMissing)
  }
  if (is.null(geno) || geno$nSNP < 2) {
    cat(paste0("Too few SNPs for window ", region, "\n"))
    return(data.frame(rho=NA, halfmax=NA, nsnp=0))
  }
  d = as.matrix(dist(geno$POS))
  d = d[upper.tri(d)]
  rsq = cor(geno$GT_minor, use="pairwise.complete.obs") ^2
  rsq = rsq[upper.tri(rsq)]

  # fits decay eqn. from Hill & Weir, via that blog post we used in brachy paper
  tryCatch({
    n = geno$nIndiv
    fit = nls(rsq ~ ((10+r*d)/((2+r*d)*(11+r*d))) *(1+((3+r*d)*(12+12*r*d+(r*d)^2))/(n*(2+r*d)*(11+r*d))),
          start=c(r=0.1), control=nls.control(maxiter=100000, warnOnly=T))
    fit = summary(fit)
    rho = fit$parameters[1]
    if (!fit$convInfo$isConv || rho < 0) {
      stop("LD decay model failed to converge (or found rho of ", rho, " < 0)")
    }
    halfmax = rho2halfmax(rho, n)
    data.frame(rho=rho, halfmax=halfmax, nsnp = geno$nSNP, stringsAsFactors = F)
  }, error = function(e) {
    cat(paste0("Got error '", as.character(e), "' in window '", region, "'\n"))
    data.frame(rho=NA, halfmax=NA, nsnp = geno$nSNP, stringsAsFactors = F)
  })
}

#' windowed_halfmax
#'
#' @param bcf Indexed BCF or VCF file
#' @param windowsize Size of each window (required unless `windows` provided)
#' @param slide Number of bases to skip between each window start (default `windowsize`)
#' @param minMAF Minimum SNP minor allle freq
#' @param maxMissing Maximum SNP missing data rate
#' @param windows restrict analyses to these windows (expected to be created with bcf_getWindows)
#' @param errors One of "remove" or "stop", see ?foreach::foreach 's .errorhandling parameter
#'
#' @return data.frame  columns rho, halfmax, region, nsnp
#' @export
windowed_halfmax = function (bcf, windowsize=NULL, slide=windowsize, minMAF=0., maxMissing=1., windows=NULL, errors="stop") {
  if (is.null(windows)) {
    if (is.null(windowsize) || is.null(slide)) stop("Invalid or missing windowsize/slide")
    cat("Making windows\n")
    windows = bcf_getWindows(bcf, windowsize=windowsize, slide = slide)
    cat(paste("Using", nrow(windows), "generated windows\n"))
  } else {
    cat(paste("Using", nrow(windows), "user-supplied windows\n"))
  }
  window_i = seq_len(nrow(windows))
  chunks = split(window_i, ceiling(window_i / (length(window_i)/1000)))
  export=c("bcf", "windows", "minMAF", "maxMissing")
  pkgs = c("dplyr", "tidyr", "purrr", "magrittr", "boringLD")
  halfmax = foreach(chunk=chunks, .combine=dplyr::bind_rows, .export = export,
                    .packages=pkgs, .errorhandling=errors) %dopar% {
    ld = purrr::map_dfr(windows[chunk,]$region, function (reg) {
            tryCatch({
              boringLD:::window_halfmax(bcf, reg, minMAF, maxMissing)
            }, error=function (err) {
              data.frame(rho=NA, halfmax=NA, nsnp=NA)
            })
          })
    dplyr::bind_cols(windows[chunk,], ld)
  }
  halfmax
}
