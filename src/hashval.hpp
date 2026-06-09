/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
***********************************************/

#pragma once

#include <cstdint>
#include <cstddef>
#include "MurmurHash3.h"

namespace GanakInt {

struct HashVal {
  uint64_t hash = 0;   // low word  -- ALSO used to index the cache hash table
  uint64_t hash2 = 0;  // high word -- extra collision discriminator
  bool operator==(const HashVal& o) const { return hash == o.hash && hash2 == o.hash2; }
  bool operator!=(const HashVal& o) const { return !(*this == o); }
};

// Returns the 128-bit MurmurHash3_x64_128 result as a HashVal.
inline HashVal murmur3_128(const void* key, const size_t len, const uint32_t seed) {
  uint64_t out[2];
  MurmurHash3_x64_128(key, (int)len, seed, out);
  return HashVal{out[0], out[1]};
}

}
