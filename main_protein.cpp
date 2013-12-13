#include <algorithm>
#include <functional>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>

#include "IndexProtein.h"
#include "IndexLocation.h"

using namespace std;

char *bwa_pg;

int main(void)
{
  IndexProtein index;
  index.process("data/protein.fasta", 1);

  char Q[] = "STD";
  char convertedQuery[123];

  std::string dna_alpha = "AGTC";
  std::string protein_alpha = "APBQCRDSETFUGVHWIYKZLXM";
  
  unsigned char nst_nt4_table[256] = {
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
    4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
  };

  int out_len = 0;
  for (int i = 0, len = strlen(Q); i < len; ++i) {
    unsigned int idx = protein_alpha.find(Q[i]);
    convertedQuery[out_len++] = dna_alpha[idx%4];
    convertedQuery[out_len++] = dna_alpha[idx/4%4];
    convertedQuery[out_len++] = dna_alpha[idx/4/4%4];
  }
  convertedQuery[out_len] = 0;
  puts(convertedQuery);
  for (int i = 0; i < out_len; ++i) {
    convertedQuery[i] = nst_nt4_table[(int)convertedQuery[i]];
  }

  printf("out_len = %d\n", out_len);

  unsigned long long numHits;
  IndexLocationList* lista = index.convertedFind(convertedQuery, out_len, /* thread id */ 0, &numHits, /* max num hits */ 10);

  printf("number of hits = %llu\n", numHits);
  printf("list size = %lu\n", lista->size());
  IndexLocationList& L = *lista;
  for (unsigned int i = 0; i < L.size(); ++i) {
    printf("seqId = %llu start = %llu rev = %d\n", L[i].sequenceId, L[i].start, (int)L[i].isReverse);
  }

  delete lista;
  return 0;
}
