#ifndef BWTINDEX_H_
#define BWTINDEX_H_

#ifdef __cplusplus
extern "C" {
#endif

  int bwa_index(int argc, char *argv[]);
  bwt_t bwt_pac2bwt(const char *fn_pac, int use_is);
  int is_bwt(ubyte_t *T, int n);

#ifdef __cplusplus
}
#endif

#endif /* BWTINDEX_H_ */

