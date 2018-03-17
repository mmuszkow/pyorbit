#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "xoroshiro128plus.h"

#define CHI_MAX  128000
#define CHI_AVG  1000.0
#define CHI_ITER (CHI_MAX * CHI_AVG)

#define POW2(x) ((x)*(x))

int chi_freq[CHI_MAX];

double chi_square(void) {
    int i;
    double chi_square = 0;
    for(i=0; i<CHI_MAX; i++) chi_square += POW2(chi_freq[i] - CHI_AVG);
    chi_square /= (double) CHI_MAX;
    return fabs(chi_square - CHI_AVG);
}

#define chi_ok(x) ((x) <= 2 * sqrt(CHI_AVG))

#define elapsed() ((clock()-start)/(float) CLOCKS_PER_SEC)

#define prng_test_start() \
    memset(chi_freq, 0, sizeof(chi_freq)); \
    start = clock();

#define prng_test_end(rand_next) \
    printf("speed=%.2fM/s\t", CHI_ITER / elapsed() / 1000000.0f); \
    chi = chi_square(); \
    printf("chi=%lf\tok=%d\t" #rand_next "\n", chi, chi_ok(chi)); \

#define prng_test(rand_next) \
    prng_test_start(); \
    for(i=0; i<CHI_ITER; i++) chi_freq[rand_next % CHI_MAX]++; \
    prng_test_end(rand_next);

int main(void) {
    int i;
    time_t start;
    double chi;
    uint64_t simd_v[4];
    xrshr128p_state_t      state;
    xrshr128p_simd_state_t simd_state;

    srand((unsigned int) time(NULL));
    xrshr128p_init(time(NULL), &state);
    xrshr128p_simd_init(time(NULL), &simd_state);
    
    prng_test(rand());
    
    prng_test(xrshr128p_next(&state));
    
    prng_test_start();
    for(i=0; i<CHI_ITER; i+=4) { /* loop is different than for other tests */
        _mm256_storeu_si256((__m256i *) simd_v, xrshr128p_simd_next(&simd_state));
        chi_freq[simd_v[0] % CHI_MAX]++;
        chi_freq[simd_v[1] % CHI_MAX]++;
        chi_freq[simd_v[2] % CHI_MAX]++;
        chi_freq[simd_v[3] % CHI_MAX]++;
    }
    prng_test_end(xrshr128p_simd_next(&simd_state));
    
    return 0;
}

