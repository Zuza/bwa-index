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
  index.process("nt.part", 1); // It will read in the index if already exists

  char Q[] = "A";
  unsigned long long numHits;
  IndexLocationList* lista = index.find(Q, strlen(Q), /* thread id */ 0, &numHits, /* max num hits */ 10);

  printf("number of hits = %llu\n", numHits);
  printf("list size = %lu\n", lista->size());
  IndexLocationList& L = *lista;
  for (unsigned int i = 0; i < L.size(); ++i) {
    printf("%llu %llu %d\n", L[i].sequenceId, L[i].start, (int)L[i].isReverse);
  }

  delete lista;
  return 0;
}
