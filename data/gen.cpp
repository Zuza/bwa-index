#include <algorithm>

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
#include <vector>
#include <iostream>

using namespace std;

#define REP(i, n) for (int i = 0; i < (n); ++i)
#define TRACE(x) cout << #x << " = " << (x) << endl

typedef long long llint;

int main(void)
{
  int n = 123123123;
  FILE* out = fopen("gen.fasta", "w");
  fputs(">dummy\n", out);
  srand(time(0));
  while (n > 0) {
    static char line[123];
    static char alpha[] = "AGTC";
    for (int i = 0; i < 80; ++i)
      line[i] = alpha[rand() % 4];
    fprintf(out, "%s\n", line);
    n -= 80;
  }
  fclose(out);
  return 0;
}
