/* Based on clhash by Daniel Lemire,Owen Kaser
 * For the original code, see: https://github.com/lemire/clhash
 *
 * LICENSE: Apache License, Version 2.0, January 2004
 *          See: http://www.apache.org/licenses/
 */

#include "stdint.h"
#include "minim.h"
#include <stdlib.h>

struct xorshift128plus_key_s {
    uint64_t part1;
    uint64_t part2;
};

typedef struct xorshift128plus_key_s xorshift128plus_key_t;

uint64_t xorshift128plus(xorshift128plus_key_t * key) {
    uint64_t s1 = key->part1;
    const uint64_t s0 = key->part2;
    key->part1 = s0;
    s1 ^= s1 << 23; // a
    key->part2 = s1 ^ s0 ^ (s1 >> 18) ^ (s0 >> 5); // b, c
    return key->part2 + s0;
}

static inline void xorshift128plus_init(uint64_t key1, uint64_t key2, xorshift128plus_key_t *key) {
    key->part1 = key1;
    key->part2 = key2;
}

void * get_random_key_for_clhash(uint64_t seed1, uint64_t seed2) {
    xorshift128plus_key_t k;
    xorshift128plus_init(seed1, seed2, &k);
    void * answer= malloc(RANDOM_BYTES_NEEDED_FOR_CLHASH);
    if (answer == 0) return 0;
    uint64_t * a64 = (uint64_t *) answer;
    for(uint32_t i = 0; i < RANDOM_64BITWORDS_NEEDED_FOR_CLHASH; ++i) {
        a64[i] =  xorshift128plus(&k);
    }
    while((a64[128]==0) && (a64[129]==1)) {
        a64[128] =  xorshift128plus(&k);
        a64[129] =  xorshift128plus(&k);
    }
    return answer;
}
