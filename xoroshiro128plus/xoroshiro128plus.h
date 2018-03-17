/* xorshiro128plus implementation (+SIMD) */

#ifndef __XOROSHIRO128PLUS_H__
#define __XOROSHIRO128PLUS_H__

#include <stdint.h>

/* this is used for expanding 64-bit seed */
uint64_t splitmix64(uint64_t x) {
	uint64_t z = (x += UINT64_C(0x9E3779B97F4A7C15));
	z = (z ^ (z >> 30)) * UINT64_C(0xBF58476D1CE4E5B9);
	z = (z ^ (z >> 27)) * UINT64_C(0x94D049BB133111EB);
	return z ^ (z >> 31);
}

typedef struct xrshr128p_state {
    uint64_t s[2];
} xrshr128p_state_t;

void xrshr128p_init(uint64_t seed, xrshr128p_state_t* state) {
    state->s[0] = seed;
    state->s[1] = splitmix64(seed);
}

#define _xrshr128p_rotl(x, k) (((x) << (k)) | ((x) >> (64 - (k))))

uint64_t xrshr128p_next(xrshr128p_state_t* state) {
	const uint64_t s0 = state->s[0];
	uint64_t s1 = state->s[1];
	const uint64_t result = s0 + s1;
	s1 ^= s0;
	state->s[0] = _xrshr128p_rotl(s0, 55) ^ s1 ^ (s1 << 14);
	state->s[1] = _xrshr128p_rotl(s1, 36);
	return result;
}

/* double in range [0,1) */
inline double xrshr128p_next_double(xrshr128p_state_t* state) {
    const uint64_t x = xrshr128p_next(state);
    const union { uint64_t i; double d; } u = { .i = UINT64_C(0x3FF) << 52 | x >> 12 };
    return u.d - 1.0;
}

#ifdef __AVX2__
#include <immintrin.h>

typedef struct xrshr128p_simd_state {
    __m256i s[2];
} xrshr128p_simd_state_t;

void xrshr128p_simd_init(uint64_t seed, xrshr128p_simd_state_t* state) {
    int i;
    uint64_t s[8];
    s[0] = seed;
    for(i=1; i<8; i++) s[i] = splitmix64(s[i-1]);
    state->s[0] = _mm256_loadu_si256((const __m256i *) s);
	state->s[1] = _mm256_loadu_si256((const __m256i *) &s[4]);
}

#define _xrshr512p_rotl(res, x, k) \
    rl1 = _mm256_slli_epi64((x), (k)); \
    rl2 = _mm256_srli_epi64((x), 64-(k)); \
    res = _mm256_or_si256(rl1, rl2);

__m256i xrshr128p_simd_next(xrshr128p_simd_state_t* state) {
    __m256i s1 = state->s[0];
    __m256i s0 = state->s[1];
    const __m256i result = _mm256_add_epi64(s0, s1);
    __m256i rot, xor_, sl, rl1, rl2;

    _xrshr512p_rotl(rot, s0, 55); 
    s1 = _mm256_xor_si256(s1, s0);
    xor_ = _mm256_xor_si256(rot, s1);
    sl = _mm256_slli_epi64(s1, 14);
    s0 = _mm256_xor_si256(xor_, sl);
    _xrshr512p_rotl(s1, s1, 36);
    
    state->s[0] = s0;
    state->s[1] = s1;
    return result;
}

#endif
#endif

