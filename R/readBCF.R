#' Reads gentoype (0, 1, 2, NA coded) and allele depths from BCF. Optionally,
#' reads only from a specific region of the reference.
#'
#' @param path Path to VCF or BCF file
#' @param region A region (like 'Chr1:1-1000') to extract. `path` must be an
#'        indexed BCF if this option is used.
#' @param samples Extract gentypes of only `samples`. Default NULL = all.
#' @param rowsAreSamples If TRUE, transpose GT and AD matricies so rows are
#'        samples and SNPs are columns. This is the opposite of a BCF.
bcf_getGTandAD = function(path, region=NULL, samples=NULL, rowsAreSamples=T, minMAF=0.0, maxMissing=1) {
  if (is.null(region)) region = ""
  if (is.null(samples)) samples = character(0)
  ret = readBCFQuery_(path, region, samples)
  nsnp = length(ret$POS)
  if (nsnp < 1) {
    warning("No SNPs in region")
    return(NULL)
  }
  # matrixify the matrices
  ret$GT = matrix(unlist(ret$GT, recursive = F), nrow=nsnp, byrow = T)
  ret$AD_ref = matrix(unlist(ret$AD_ref, recursive = F), nrow=nsnp, byrow = T)
  ret$AD_alt = matrix(unlist(ret$AD_alt, recursive = F), nrow=nsnp, byrow = T)
  # Adjust missing gentoypes: -1 is missing in C++ output
  ret$GT[ret$GT<0] = NA
  # Missing rate and AF
  N = ncol(ret$GT)
  ret$nIndiv = N
  ret$nSNP = nsnp
  ret$MissRate = rowSums(is.na(ret$GT))/N
  ret$AF = rowSums(ret$GT, na.rm = T)/(2*rowSums(!is.na(ret$GT)))
  ret$AF[!is.finite(ret$AF)] = 0 # as one can get NANs due to div by 0!
  # Recode 012 to be minor allele counts
  ret$GT_minor = ret$GT
  ret$GT_minor[ret$AF>0.5,] = 2 - ret$GT_minor[ret$AF>0.5,]
  ret$MAF = rowSums(ret$GT_minor, na.rm = T)/(2*rowSums(!is.na(ret$GT_minor)))
  ret$Ho = rowSums(ret$GT_minor==1, na.rm = T)/(2*rowSums(!is.na(ret$GT_minor)))
  # Postions are zero-based in BCF, we want 1-based
  ret$POS = ret$POS + 1
  snp.keep = ret$MAF >= minMAF & ret$MissRate <= maxMissing
  ret$nSNP = sum(snp.keep)
  for (mat in c("GT", "AD_ref", "AD_alt", "GT_minor")) {
    # Keep good SNPs
    ret[[mat]] = ret[[mat]][snp.keep,]
    colnames(ret[[mat]]) = ret$Samples
    if (rowsAreSamples) {
      ret[[mat]] = t(ret[[mat]])
    }
  }
  for (vec in c("CHROM", "POS", "AF", "MissRate", "MAF", "Ho")) {
    ret[[vec]] = ret[[vec]][snp.keep]
  }
  ret
}


#' Reads contig names and lengths from BCF file as a data frame.
#'
#' @param path Path to VCF or BCF file
bcf_getContigs = function(path) {
  ret = readBCFContigs_(path)
  # Adjust missing gentoypes: -1 is missing in C++ output
  as.data.frame(ret, stringsAsFactors=F)
}


#' Generate a list of regions which slide over the genome
#'
#' @param path Path to VCF or BCF file
#' @param windowsize Size of each window, in bases
#' @param slide Step size of each window, in bases
#' @param chrom Only generate windows on contig 'chrom'
#' @param from Only generate windows starting at `from` on `chrom`
#' @param to Only generate windows until `to` on `chrom`
#'
#' @export
bcf_getWindows = function(path, windowsize=1000, slide=windowsize, chrom=NULL, from=NULL, to=NULL) {
  contigs = bcf_getContigs(path)
  res = NULL
  if (is.null(chrom) && !is.null(from) && !is.null(to)) {
    stop("chrom must be specified if from and to are specified")
  }
  for (ctg in seq_len(nrow(contigs))) {
    name = contigs$names[ctg]
    if (!is.null(chrom) && name != chrom) {
      next
    }
    ctglen = contigs$lengths[ctg]
    if (!is.null(chrom) && !is.null(from) && !is.null(to)) {
      from = min(from, ctglen)
      to = min(to, ctglen) - windowsize + slide
      if (from >= to) {
        stop("from must be < to-windowsize")
      }
    } else {
      from = 1
      to = ctglen - windowsize + slide
    }
    for (start in seq(from, to, slide)) {
      stop = start + windowsize - 1
      region = sprintf("%s:%d-%d", name, start, stop)
      res = rbind(res, data.frame(region,contig=name, start, stop, stringsAsFactors = F))
    }
  }
  res
}
