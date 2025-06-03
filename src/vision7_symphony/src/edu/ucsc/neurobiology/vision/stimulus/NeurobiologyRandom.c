#include <time.h>

// This replicates the original Macintosh Toolbox Random() routine.
// The funciton is defined here explicitly for portability
// and independence from changes in the MacOS.  EJC 1999-12-22
short StimRandShort(int *seed) {

  int temp1, temp2, temp3;
  short result;
  
  temp1 = (*seed & 0xFFFF) * 0x41A7;
  temp2 = (*seed >> 16) * 0x41A7 + (temp1 >> 16);
  temp3 = temp2 * 2 >> 16;
  temp1 = (temp1 & 0xFFFF) - 0x7FFFFFFF;
  temp2 &= 0x7FFF;
  temp1 += (temp2 << 16) | (temp2 >> 16) + temp3;
  if (temp1 < 0) temp1 += 0x7FFFFFFF;
  *seed = temp1;

  result = temp1 & 0xFFFF;
  if (((unsigned short)result) == 0x8000) result = 0;
  return result;
}

int main(int argc, char *argv[]) {

  int i = 0;
  int seed = 0;
  int limit = 0;
  short value = 0;
  clock_t start;
  clock_t end;
  double diff;

  sscanf(argv[1],"%i",&seed);
  sscanf(argv[2],"%i",&limit);
  printf("%i %i\n",seed,limit);

  start = clock();
  for (i=0; i<limit; i++) {
    value = StimRandShort(&seed);
  }
  end = clock();
  printf("%i %i %i\n",seed,limit,value);

  diff = (end-start)*1000./CLOCKS_PER_SEC;
  printf("time: %f\n",diff);
}

