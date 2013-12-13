#include <algorithm>
#include <functional>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>

#include "IndexBWA.h"
#include "IndexLocation.h"

using namespace std;

char *bwa_pg;

int main(void)
{
  IndexBWA index;
  //  index.process("data/ls_orchid.fasta", 1); // It will read in the index if already exists
  index.process("data/gen.fasta", 1);

  char Q[] = "AGTCAAAGTCCGAGATAA";
  char convertedQuery[123];

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

  unsigned long long suma = 0;

  for (int i = 0; i < 1000000; ++i) {
    const int len = 20;
    static char alphabet[] = "AGTC";
    for (int i = 0; i < len; ++i)
      convertedQuery[i] = nst_nt4_table[(int)alphabet[rand() % 4]];
    unsigned long long numHits;
    IndexLocationList* lista = index.convertedFind(convertedQuery, strlen(Q), /* thread id */ 0, &numHits, /* max num hits */ 10);
    IndexLocationList& L = *lista;
    // for (unsigned int i = 0; i < L.size(); ++i) {
    //   printf("seqId = %llu start = %llu rev = %d\n", L[i].sequenceId, L[i].start, (int)L[i].isReverse);
    // }
    suma += numHits;
    delete lista;
  }

  printf("number of hits = %llu\n", suma);
  return 0;
}
