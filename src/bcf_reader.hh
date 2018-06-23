#include <exception>
#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <htslib/hts.h>
#include <htslib/vcf.h>

namespace KDMBCF {

using namespace std;

class BCFReader
{
public:
  template <typename T> using Mat2d = vector<vector<T>>;
  BCFReader (const string &filename_) {
    filename = filename_;
    if((bcf = bcf_open(filename.c_str(), "r")) == NULL) {
      throw runtime_error("Unable to open file.");
    }
    if((header = bcf_hdr_read(bcf)) == NULL) {
      throw runtime_error("Unable to read header.");
    }
    record = bcf_init();
    nsamp = bcf_hdr_nsamples(header);
  }
  ~BCFReader ()
  {
    bcf_hdr_destroy(header);
    bcf_destroy(record);
    bcf_close(bcf);
  }

  vector<string> get_sample_names()
  {
      vector<string> samples;
      for (int32_t i = 0; i < nsamp; i++) {
          samples.push_back(header->samples[i]);
      }
      return samples;
  }

  bool set_samples(const vector<string> &samples)
  {
      string samplecsv = "";
      for (size_t i = 0; i < samples.size(); i++) {
          samplecsv += samples[i];
          if (i < samples.size() - 1)
              samplecsv += ",";
      }
      if (bcf_hdr_set_samples(header, samplecsv.c_str(), 0) != 0) {
          return false;
      }
      nsamp = bcf_hdr_nsamples(header);
      if (nsamp != samples.size()) return false;
      return true;
  }

  bool get_contig_names_lengths(vector<string> &names, vector<int32_t> &lengths)
  {
    int32_t nctg = header->n[BCF_DT_CTG];
    names.clear();
    lengths.clear();
    for (int32_t i = 0; i < nctg; i++) {
      bcf_idpair_t *ctg = header->id[BCF_DT_CTG];
      names.emplace_back(ctg[i].key);
      lengths.emplace_back(ctg[i].val->info[0]);
    }
    return true;
  }

  bool read_region(const string &region)
  {
    CHROM.clear();
    POS.clear();
    GT.clear();
    AD_ref.clear();
    AD_alt.clear();
    hts_idx_t *idx = bcf_index_load(filename.c_str());
    hts_itr_t *iter;
    iter = bcf_itr_querys(idx, header, region.c_str());
    while(bcf_itr_next(bcf, iter, record) == 0) {
      process_record();
    }
    bcf_itr_destroy(iter);
    hts_idx_destroy(idx);
  }

  bool read_all()
  {
    CHROM.clear();
    POS.clear();
    GT.clear();
    AD_ref.clear();
    AD_alt.clear();
    while(bcf_read(bcf, header, record) == 0) {
      process_record();
    }
  }

  vector< string > CHROM;
  vector< uint64_t > POS;
  Mat2d<int32_t> GT;
  Mat2d<int32_t> AD_ref;
  Mat2d<int32_t> AD_alt;

protected:
  htsFile *bcf = NULL;
  bcf_hdr_t *header = NULL;
  bcf1_t *record = NULL;
  int32_t *buffer = NULL;
  int32_t buffersz = 0;
  int32_t nsamp = 0;
  string filename;

  void process_record()
  {
      vector<int32_t> GT1, AD_ref1, AD_alt1;
      if (!bcf_is_snp(record) ||   record->n_allele > 2) return;
      if (!get_GT(GT1)) return;
      if (!get_AD(AD_ref1, AD_alt1)) return;
      GT.push_back(GT1);
      AD_ref.push_back(AD_ref1);
      AD_alt.push_back(AD_alt1);
      CHROM.emplace_back(bcf_hdr_id2name(header, record->rid));
      POS.push_back(record->pos);
  }

  int32_t fill_buffer(const char *tag) {
    int32_t nentry = bcf_get_format_int32(header, record, tag, &buffer, &buffersz);
    int32_t ploid = nentry/nsamp;
    if (nentry < 0) return nentry;
    if (ploid != 2) return 0;
    return nentry;
  }

  bool get_GT(vector<int32_t> &GT1)
  {
    int32_t n = fill_buffer("GT");
    if (n < 1) return false;
    GT1.resize(nsamp);
    for (int32_t i=0; i<nsamp; i++) {
      int32_t a=bcf_gt_allele(buffer[i]), b=bcf_gt_allele(buffer[i+1]);
      int32_t altcnt = (a<0 || b< 0) ? -1 : a + b;
      GT1[i] = altcnt;
    }
    return true;
  }

  bool get_AD(vector<int32_t> &AD_ref1, vector<int32_t> &AD_alt1)
  {
    int32_t n = fill_buffer("AD");
    int32_t ploid = n/nsamp;
    if (n < 1) return false;

    AD_ref1.resize(nsamp);
    AD_alt1.resize(nsamp);
    for (int32_t i=0; i<nsamp; i++) {
      int32_t *ADs = buffer + i*ploid;
      AD_ref1[i] = ADs[0];
      AD_alt1[i] = ADs[1];
    }
    return true;
  }

};


} /* namespace gusld */
