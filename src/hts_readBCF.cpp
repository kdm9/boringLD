#include <Rcpp.h>
#include <string>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <vector>
#undef KDM_BCF_READER_USE_AD
#include "bcf_reader.hh"

using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins(cpp11)]]


// [[Rcpp::export]]
List readBCFQuery_(SEXP fname, SEXP reg, SEXP samplenames) {
  string filename = as<string>(fname);
  string region = as<string>(reg);
  vector<string> request_samples = as<vector<string>>(samplenames);
  KDMBCF::BCFReader rdr(filename);

  if (request_samples.size() > 0) {
      if (!rdr.set_samples(request_samples)) return NULL;
  }
  vector<string> samples = rdr.get_sample_names();
  if (region == "") {
    if (!rdr.read_all()) return NULL;
  } else {
    if (!rdr.read_region(region)) return NULL;
  }

  return List::create(
    Named("CHROM")=wrap(rdr.CHROM),
    Named("POS")=wrap(rdr.POS),
    Named("GT")=wrap(rdr.GT),
    Named("Samples")=wrap(samples)
#ifdef KDM_BCF_READER_USE_AD
    , Named("AD_ref")=wrap(rdr.AD_ref)
    , Named("AD_alt")=wrap(rdr.AD_alt)
#endif
    );
}

// [[Rcpp::export]]
List readBCFContigs_(SEXP fname) {
  string filename = as<string>(fname);
  KDMBCF::BCFReader rdr(filename);
  vector<string> ctg_names;
  vector<int32_t> ctg_lenghts;

  rdr.get_contig_names_lengths(ctg_names, ctg_lenghts);

  return List::create(Named("names")=wrap(ctg_names),
                      Named("lengths")=wrap(ctg_lenghts));
}

// [[Rcpp::export]]
CharacterVector readBCFSamples_(SEXP fname) {
  KDMBCF::BCFReader rdr(as<string>(fname));

  return wrap(rdr.get_sample_names());
}
