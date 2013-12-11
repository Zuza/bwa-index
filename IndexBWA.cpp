/*
 * IndexBWA.cpp
 *
 *  Created on: Oct 13, 2013
 *      Author: ivan
 */

#include "IndexBWA.h"
#include <iostream>


IndexBWA::IndexBWA()
{
	idx_ = NULL;
	itrs_ = NULL;
	numItrs_ = 0;
}

IndexBWA::~IndexBWA()
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

// IndexBWA::IndexBWA(SequenceSet *sequences, unsigned int seedLength)
// {
// 	ErrorHandling::functionNotImplemented(((std::string) __FUNCTION__));
// }

// IndexBWA::IndexBWA(std::string sequencesPath)
// {
// 	ErrorHandling::functionNotImplemented(((std::string) __FUNCTION__));
// }

void IndexBWA::process(std::string sequencesPath, unsigned int numItrs)
{
	char arg0[] = "index";
	char *arg1 = (char *) sequencesPath.c_str();

    char* argv[] = { &arg0[0], &arg1[0], NULL };
    int   argc   = (int)(sizeof(argv) / sizeof(argv[0])) - 1;

    FILE *fp=NULL;
    fp = fopen((sequencesPath + ((std::string) ".bwt")).c_str(), "r");
    if (fp == NULL)
    	bwa_index(argc, argv);
    else {
    	std::cout << "BWT index already exists, skipping index creation. Loading index from file." << std::endl;
      fclose(fp);
    }

	idx_ = bwa_idx_load(sequencesPath.c_str(), BWA_IDX_ALL); // load the BWA index

	numItrs_ = numItrs;
	itrs_ = (smem_i **) calloc(numItrs_, sizeof(smem_i *));
	for (unsigned int i=0; i<numItrs_; i++)
		itrs_[i] = smem_itr_init(idx_->bwt);
}

bwaidx_t *IndexBWA::getIndex()
{
	return idx_;
}

IndexLocationList* IndexBWA::find(char *query, unsigned int queryLength, unsigned int threadId, unsigned long long int *retNumHits, unsigned long int maxNumHits)
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
		for (i = 0; i < a->n; ++i)
		{
			bwtintv_t *p = &a->a[i];

			if ((uint32_t)p->info - (p->info>>32) < queryLength)
				continue;

			numHits += p->x[2];

			// Check the number of hits.
			if (p->x[2] <= maxNumHits)
			{
				ret->resize((ret->size() + p->x[2]));

				for (k=0; k < p->x[2]; ++k)
				{
					bwtint_t pos;
					int len, is_rev, ref_id;

					len  = (uint32_t)p->info - (p->info>>32);
					pos = bns_depos(idx_->bns, bwt_sa(idx_->bwt, p->x[0] + k), &is_rev);

					if (is_rev)
					{
						pos -= len - 1;
					}

					bns_cnt_ambi(idx_->bns, pos, len, &ref_id);

					// Note: the hit position (pos) is absolute, which means that all reference sequences that were indexed in the original FASTA file were first truncated into one "FASTA sequence", so the location that we get is absolute with respect to the beginning of the FASTA file. That's why we need to subtract the offset of the current reference sequence that the hit is located on.
					IndexLocation currentLocation(ref_id, (pos - idx_->bns->anns[ref_id].offset), (unsigned char) (is_rev==0));
					(*ret)[numRet] = currentLocation;

					numRet += 1;
				}
			}
		}
	}

	*retNumHits = numHits;

	return ret;
}

// void IndexBWA::verbose(std::ostream &outStream)
// {
// 	ErrorHandling::functionNotImplemented(((std::string) __FUNCTION__));
// }

std::string IndexBWA::getHeader(unsigned long long int sequenceId)
{
	return ((std::string) idx_->bns->anns[sequenceId].name);
}

smem_i** IndexBWA::getItrs()
{
	return itrs_;
}
