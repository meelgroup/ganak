/******************************************
Copyright (C) 2026 Randal Bryant
                   Minor changes by Mate Soos

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


#include <stdbool.h>
#include <gmp.h>
#include <mpfr.h>
#include <mpfi.h>

#pragma once

/* Allow this headerfile to define C++ constructs if requested */
#ifdef __cplusplus
#define CPLUSPLUS
#endif

#ifdef CPLUSPLUS
extern "C" {
#endif


typedef struct {
    unsigned prec;
    unsigned qsize;  // bytes for mpq, 0 for mpfi
    mpfi_t mval;
    mpq_t qval;
} mpqi_t, *mpqi_ptr;

unsigned mpqi_get_default_prec();
void mpqi_set_default_prec(unsigned prec);

void mpqi_reset();
void mpqi_set_parameters(size_t initial_bytes, size_t final_bytes, unsigned long cross_count);
void mpqi_get_parameters(size_t *initial_bytes, size_t *final_bytes, unsigned long *cross_count);

void mpqi_init2(mpqi_ptr mp, unsigned prec);
void mpqi_init(mpqi_ptr mp);
void mpqi_clear(mpqi_ptr mp);

void mpqi_set_d(mpqi_ptr mp, double d);
void mpqi_set_q(mpqi_ptr mp, mpq_srcptr q);
void mpqi_set_m(mpqi_ptr mp, mpfi_srcptr m);
void mpqi_set(mpqi_ptr mp, mpqi_ptr m);
double mpqi_get_d(mpqi_ptr mp);

void mpqi_left(mpfr_ptr mid, mpqi_ptr mp);
void mpqi_right(mpfr_ptr mid, mpqi_ptr mp);

void mpqi_mid(mpfr_ptr mid, mpqi_ptr mp);
void mpqi_mid_q(mpq_ptr mid, mpqi_ptr mp);



void mpqi_mul_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q);
void mpqi_mul_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr m);
void mpqi_mul(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2);

void mpqi_add_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q);
void mpqi_add_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr m);
void mpqi_add(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2);


bool mpqi_is_zero(mpqi_ptr mp);

double digit_precision_mpqi(mpqi_ptr mp);

#ifdef CPLUSPLUS
}
#endif
