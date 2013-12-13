/*
 * IndexProtein.cpp
 *
 */

#include "IndexProtein.h"
#include <iostream>

IndexProtein::IndexProtein()
{
	idx_ = NULL;
	itrs_ = NULL;
	numItrs_ = 0;
}

IndexProtein::~IndexProtein()
{
	if (idx_)
		bwa_idx_destroy(idx_);

	if (itrs_)
	{
		for (unsigned int i=0; i<numItrs_; i++)
			smem_itr_destroy(itrs_[i]);

		free(itrs_);
	}
}

void IndexProtein::process(std::string sequencesPath, unsigned int numItrs)
{
	char arg0[] = "index";
  std::string fullPath = sequencesPath + ((std::string) ".conv");
	char *arg1 = (char *) fullPath.c_str();

  char* argv[] = { &arg0[0], &arg1[0], NULL };
  int   argc   = (int)(sizeof(argv) / sizeof(argv[0])) - 1;

  std::string dna_alpha = "AGTC";
  std::string protein_alpha = "APBQCRDSETFUGVHWIYKZLXM";

  FILE *fp=NULL;
  fp = fopen(arg1, "r");
  if (fp == NULL) {
    static char line[16000];
    static char outline[16000*3];
    FILE* in = fopen(sequencesPath.c_str(), "r");

    fp = fopen(arg1, "w");
    bool comment = false;
    while (fgets(line, 16000, in)) {
      if (line[0] == '>') comment = true;
      int len = strlen(line);
      if (comment) {
        fputs(line, fp);
        if (line[len-1] == '\n') comment = false;
      } else {
        int out_len = 0;
        for (int i = 0; i < len; ++i)
          if (isspace(line[i]))
            outline[out_len++] = line[i];
          else {
            unsigned int idx = protein_alpha.find(line[i]);
            outline[out_len++] = dna_alpha[idx%4];
            outline[out_len++] = dna_alpha[idx/4%4];
            outline[out_len++] = dna_alpha[idx/4/4%4];
          }
        fputs(outline, fp);
      }
    }

    fclose(in);
    fclose(fp);

    bwa_index(argc, argv);
  }
  else {
    std::cout << "BWT index already exists, skipping index creation. Loading index from file." << std::endl;
    fclose(fp);
  }

	idx_ = bwa_idx_load(fullPath.c_str(), BWA_IDX_ALL); // load the BWA index

	numItrs_ = numItrs;
	itrs_ = (smem_i **) calloc(numItrs_, sizeof(smem_i *));
	for (unsigned int i=0; i<numItrs_; i++)
		itrs_[i] = smem_itr_init(idx_->bwt);
}

bwaidx_t *IndexProtein::getIndex()
{
	return idx_;
}

IndexLocationList* IndexProtein::convertedFind(char *query, unsigned int queryLength, unsigned int threadId, unsigned long long int *retNumHits, unsigned long int maxNumHits)
{
	unsigned long long int numHits=0;
	unsigned long long int i=0, k=0;
	unsigned long long int numRet=0;
	IndexLocationList *ret=NULL;

	smem_i *itr_=itrs_[threadId];
	ret = new IndexLocationList;

	smem_set_query(itr_, queryLength, (uint8_t*) query);

	const bwtintv_v *a;
	while ((a = smem_next(itr_, queryLength, 0)) != 0)
	{
		for (i = 0; i < a->n; ++i) // ?? what does a->n represent
		{
			bwtintv_t *p = &a->a[i];

			if ((uint32_t)p->info - (p->info>>32) < queryLength) // ??
				continue;

			// Check the number of hits.
			if (numRet <= maxNumHits)
			{
				for (k=0; k < p->x[2] && numRet <= maxNumHits; ++k)
				{
					bwtint_t pos;
					int len, is_rev, ref_id;

					len  = (uint32_t)p->info - (p->info>>32);
					pos = bns_depos(idx_->bns, bwt_sa(idx_->bwt, p->x[0] + k), &is_rev);

          // if it's reversed, skip it because it doesn't exist in protein sequences
					if (is_rev) {
            continue;
					}
					bns_cnt_ambi(idx_->bns, pos, len, &ref_id);

          // only positions 0, 3, 6, ... correspond to real protein starts
          if ((pos - idx_->bns->anns[ref_id].offset) % 3 != 0) {
            continue;
          }

					// Note: the hit position (pos) is absolute, which means that
          // all reference sequences that were indexed in the original
          // FASTA file were first truncated into one "FASTA sequence",
          // so the location that we get is absolute with respect to the 
          // beginning of the FASTA file. That's why we need to subtract
          // the offset of the current reference sequence that the hit is
          // located on.
					IndexLocation currentLocation(ref_id, (pos - idx_->bns->anns[ref_id].offset) / 3, (unsigned char)is_rev);
					(*ret).push_back(currentLocation);

          numHits += 1;
					numRet += 1;
				}
			}
		}
	}

	*retNumHits = numHits;
	return ret;
}

// void IndexProtein::verbose(std::ostream &outStream)
// {
// 	ErrorHandling::functionNotImplemented(((std::string) __FUNCTION__));
// }

std::string IndexProtein::getHeader(unsigned long long int sequenceId)
{
	return ((std::string) idx_->bns->anns[sequenceId].name);
}

smem_i** IndexProtein::getItrs()
{
	return itrs_;
}
