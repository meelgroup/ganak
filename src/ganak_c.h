/******************************************
Copyright (C) 2024 Authors of GANAK, see AUTHORS file

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

/*
 * Pure C ABI for Ganak.
 *
 * Use together with arjun_c.h: build a SimplifiedCNF, run Arjun preprocessing,
 * then hand the resulting cnf to ganak_counter_setup_from_simpcnf() and call
 * ganak_counter_count() to obtain the (projected, optionally weighted) model
 * count as an arjun_field_t.
 */

#ifndef GANAK_C_H
#define GANAK_C_H

#include <stdint.h>
#include <stddef.h>
#include <arjun/arjun_c.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Opaque type. */
typedef struct ganak_counter ganak_counter_t;

/* Returns the message from the most recent failing Ganak C call on the
 * current thread, or "" if none. Owned by ganak; do not free. */
const char* ganak_last_error(void);

/* Construct a counter with a FieldGen (same one used to build the CNF).
 *   verb:  verbosity (0 = silent).
 *   seed:  RNG seed.
 *   fg:    field generator (MPQ or MPZ); not consumed, caller still owns. */
ganak_counter_t* ganak_counter_new(uint32_t verb, uint64_t seed,
                                   const arjun_fgen_t* fg);

void ganak_counter_free(ganak_counter_t* c);

/* Mirrors setup_ganak<>(simp_cnf, counter): copy vars, indep support, optional
 * indep support, weights and clauses from the SimplifiedCNF into the counter.
 * Must be called exactly once before ganak_counter_count(). */
int ganak_counter_setup_from_simpcnf(ganak_counter_t* c,
                                     const arjun_simpcnf_t* cnf);

/* Run the counting algorithm. Returns a NEWLY-OWNED arjun_field_t with the
 * count (caller must arjun_field_free it). On error, returns NULL and sets
 * ganak_last_error().
 *
 * num_threads: 1 = single-threaded (matches example.cpp). */
arjun_field_t* ganak_counter_count(ganak_counter_t* c,
                                   uint8_t bits_threads,
                                   int num_threads,
                                   int debug_threads);

/* True if the count is approximate (e.g. ApproxMC fallback was used). */
int ganak_counter_is_approximate(const ganak_counter_t* c);

#ifdef __cplusplus
}  /* extern "C" */
#endif

#endif /* GANAK_C_H */
