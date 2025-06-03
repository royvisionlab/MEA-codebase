#include <stdio.h>
#include <math.h>

static long A = 0x5DEECE66DL;
static long B = 0xBL;
static long mask = (1L << 48) - 1;

long init(long s) {
  return (s ^ A) & mask;
}

long advance(long s) { 
  return (A*s + B) & mask;
}

long advance_direct(long s, int n) {
  long Ad = 1, Bd = 0;
  for (; n > 0; --n) {
    Ad *= A;
    Bd = Bd*A + 1;
  }
  Bd *= B;
  return (Ad*s + Bd) & mask;
}

long advance_iterate(long s, int n) {
  for (; n > 0; --n) s = advance(s);
  return s;
}

int main(void) {
  long s;
  int n = 1000;

  s = init(11111);
  s = advance_iterate(s, n);
  printf("%ld\n", s);
  
  s = init(11111);
  s = advance_direct(s, n);
  printf("%ld\n", s);
  
  return 0;
}
