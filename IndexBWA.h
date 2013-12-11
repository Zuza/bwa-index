/*
 * IndexBWA.h
 *
 *  Created on: Oct 13, 2013
 *      Author: ivan
 */

#ifndef INDEXBWA_H_
#define INDEXBWA_H_

#include <stdio.h>
#include <string.h>
#include <vector>
#include <list>
#include "bwa/bwamem.h"
#include "bwa/kseq.h" // for the FASTA/Q parser
#include "bwa/bwt.h"
#include "bwa/bwa.h"
#include "bwa/bwase.h"
#include "bwa/utils.h"

#include "bwtindex.h"
#include "IndexLocation.h"

typedef std::vector<IndexLocation> IndexLocationList;

class IndexBWA
{
public:
	IndexBWA();
  //	IndexBWA(SequenceSet *sequences, unsigned int seedLength);
  //	IndexBWA(std::string sequencesPath);
	~IndexBWA();

	/** Process a fasta file to create an index. If an index with the given filename already exists, function skips index creation.
	 * @param sequencesPath Path to an input FASTA file.
	 * @param numItrs Number of threads that will be used (paralelizirao sam svoj aligner s OpenMP-om, pa sam ove iteratore multiplicirao tako da imam za svaku jezgru zaseban. Cini mi se da je radilo dosta brze nego da sam stavio pragma omp critical, ali opet, ovo je moje prvo iskustvo s paralelizacijom. Ako koristis jednu jezgru, ovaj dio koji koristi numIters mozes maknuti (i smem_i **itrs_; dolje u private dijelu moze biti onda jednostruki pointer).
	 */
	void process(std::string sequencesPath, unsigned int numItrs);

	/** Finds all exact matches of a query (given as a char pointer and its length), and returns them as a list (or a vector, depends on the typedef at the begining of this file).
	 * @param query Char pointer to the beginning of a seed. Note that the query sequence needs to be 2-bit coded, as used by BWA, otherwise it won't work correctly.
	 * @param queryLength Length of the seed. Explicitly stated
	 * @param threadId ID of the thread that currently accesses the index (thread ID is used as the index for the itrs_ array).
	 * @param retNumHits Total number of exact matches for a given seed.
	 * @param maxNumHits Maximum allowed number of hits. Sometimes the number of hits can be really large (dozens of thousands), which can indicate highly repetitive regions, and may introduce noise in the alignment procedure. Only if the number of hits is lesser than (or equal to) maxNumHits will the hit locations be returned. Otherwise, the return value will be equal to NULL;
	 * @return Pointer to the list of hit locations. Equal to NULL if there are no hits, or the number of hits is larger than maxNumHits.
	 */
	IndexLocationList* find(char *query, unsigned int queryLength, unsigned int threadId, unsigned long long int *retNumHits, unsigned long int maxNumHits);
	/** Returns the header of an indexed reference sequence.
	 * @param sequenceId ID of the sequence.
	 * @return String containing the header of the sequence.
	 */
	std::string getHeader(unsigned long long int sequenceId);

  //	void verbose(std::ostream &outStream);

	bwaidx_t *getIndex();
	smem_i** getItrs();

private:
	bwaidx_t *idx_;
	smem_i **itrs_;
	unsigned int numItrs_;
};

#endif /* INDEXBWA_H_ */
