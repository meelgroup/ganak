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
#include <vector>

namespace GanakInt {

// Weisfeiler-Lehman canonical component information.
// Computed once in CompAnalyzer for small components (nVars <= wl_canonize_threshold).
// Invariant to variable renaming: two structurally isomorphic components produce identical
// sorted_canon_clauses (and therefore the same hash), enabling cache hits across
// different variable numberings.
struct CanonInfo {
  bool valid = false;

  // canon_vars[i] = original var_id of the i-th canonical variable.
  // Variables are sorted by WL color (degree-based initial coloring + 1 WL round).
  // Tiebreaker: original var_id for stability.
  std::vector<uint32_t> canon_vars;

  // sorted_canon_clauses[j] = sorted list of canonical var indices (0..nVars-1) in clause j.
  // The outer vector is sorted lexicographically so the representation is canonical.
  // Binary clauses appear as 2-element inner vectors.
  std::vector<std::vector<uint32_t>> sorted_canon_clauses;

  // Pre-computed structural hash derived from sorted_canon_clauses and nVars.
  // Used directly as the cache key hash for both HashedComp and DiffPackedComp modes.
  uint64_t hash = 0;
};

} // namespace GanakInt
