/* Based on clhash by Daniel Lemire,Owen Kaser
 * For the original code, see: https://github.com/lemire/clhash
 *
 * LICENSE: Apache License, Version 2.0, January 2004
 *          See: http://www.apache.org/licenses/
 */
#ifndef MINIM_CLHASH__
#define MINIM_CLHASH__

#include "stdint.h"

#ifdef __cplusplus
extern "C" {
#endif

enum {RANDOM_64BITWORDS_NEEDED_FOR_CLHASH=133,RANDOM_BYTES_NEEDED_FOR_CLHASH=133*8};
void * get_random_key_for_clhash(uint64_t seed1, uint64_t seed2);

#ifdef __cplusplus
} // extern "C"
#endif
#endif
