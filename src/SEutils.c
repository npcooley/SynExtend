#include "SEutils.h"

int *sample(int n){
  int *r = malloc(sizeof(int) * n);
  int j;
  for (int i=0; i<n; i++){
    j = irand() % (i+1);
    if (j != i)
      r[i] = r[j];
    r[j] = i;
  }

  return r;
}

void shuffle_int_(int *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_uint_(uint *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}


void shuffle_double_(double *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void shuffle_char_(char *x, int n){
  int j, tmp;
  for (int i=(n-1); i>0; i--){
    j = irand() % (i+1);
    tmp = x[j];
    x[j] = x[i];
    x[i] = tmp;
  }
}

void seedRNGState64(struct RNGstate64 *r, uint64_t seed){
  // splitmix64 initialization, see https://en.wikipedia.org/wiki/Xorshift
  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[0] = (uint32_t) seed;
  r->state[1] = (uint32_t) (seed >> 32);
  return;
}

uint64_t xorshift128p(struct RNGstate64 *r){
  uint64_t t = r->state[0];
  uint64_t const s = r->state[1];

  r->state[0] = s;
  t ^= (t << 23);
  t ^= (t >> 18);
  t ^= s ^ (s >> 5);
  r->state[1] = t;

  return t+s;
}

void seedRNGState32(struct RNGstate32 *r, uint64_t seed){
  // splitmix64 initialization, see https://en.wikipedia.org/wiki/Xorshift
  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[0] = (uint32_t) seed;
  r->state[1] = (uint32_t) (seed >> 32);

  seed += 0x9E3779B97F4A7C15;
  seed = (seed ^ (seed >> 30)) * 0xBF58476D1CE4E5B9;
  seed = (seed ^ (seed >> 27)) * 0x94D049BB133111EB;
  seed ^= (seed >> 31);
  r->state[2] = (uint32_t) seed;
  r->state[3] = (uint32_t) (seed >> 32);
  return;
}

// adapted from https://prng.di.unimi.it/xoshiro128plus.c
uint32_t xorshift32b(struct RNGstate32 *rngstate) {
  uint32_t *r = rngstate->state;
  const uint32_t result = r[0] + r[3];
  const uint32_t t = r[1] << 9;

  r[2] ^= r[0];
  r[3] ^= r[1];
  r[1] ^= r[2];
  r[0] ^= r[3];

  r[2] ^= t;

  r[3] = rotl(r[3], 11);

  return result;
}


/* needs C11 :/
#define shuffle(x, n) _Generic((x),
  int *: shuffle_i(x, n),
  double *: shuffle_d(x, n),
  char *: shuffle_c(x, n))
*/