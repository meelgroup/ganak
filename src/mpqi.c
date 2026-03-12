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

#include "mpqi.h"

#include <stdio.h>
#include <limits.h>
#include <stdatomic.h>
#include <mpfi.h>

#define MAX_DIGIT_PRECISION 1e6

static unsigned mpqi_default_prec = 64;

/* Limits on byte size for MPQ */
static int initial_max_size = 10000;
static int final_max_size = 1000;
/* When to transition */
static unsigned long crossover_opcount = 1000000;

/* Current status */
static _Atomic unsigned long opcount = 0;

unsigned mpqi_get_default_prec() {
    return mpqi_default_prec;
}

void mpqi_set_default_prec(unsigned prec) {
    mpqi_default_prec = prec;
}

void mpqi_reset() {
    opcount = 0;
}

void mpqi_set_parameters(size_t initial_bytes, size_t final_bytes, unsigned long cross_count) {
    initial_max_size = initial_bytes;
    final_max_size = final_bytes;
    crossover_opcount = cross_count;
}

void mpqi_get_parameters(size_t *initial_bytes, size_t *final_bytes, unsigned long *cross_count) {
    *initial_bytes = initial_max_size;
    *final_bytes = final_max_size;
    *cross_count = crossover_opcount;
}


static size_t mpq_bytes(mpq_srcptr val) {
    /* Overhead */
    size_t size = 32;
    mpz_t mp_num, mp_den;
    mpz_init(mp_num); mpz_init(mp_den);
    mpq_get_num(mp_num, val);
    size += mpz_size(mp_num) * sizeof(mp_limb_t);
    mpq_get_den(mp_den, val);
    size += mpz_size(mp_den) * sizeof(mp_limb_t);
    mpz_clear(mp_num); mpz_clear(mp_den);
    return size;
}

void mpqi_init2(mpqi_ptr mp, unsigned prec) {
    mp->prec = prec;
    mpq_init(mp->qval);
    mpq_set_d(mp->qval, 0.0);
    mp->qsize = mpq_bytes(mp->qval);
}

void mpqi_init(mpqi_ptr mp) {
    mpqi_init2(mp, mpqi_default_prec);
}

void mpqi_clear(mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpq_clear(mp->qval);
    else
        mpfi_clear(mp->mval);
}

static bool size_ok(unsigned bytes) {
    unsigned max_size = (unsigned)
        (opcount < crossover_opcount ? initial_max_size : final_max_size);
    return bytes <= max_size;
}

static void mpqi_arg_check(mpqi_ptr mp) {
    if (size_ok(mp->qsize))
        return;
    mpfi_init2(mp->mval, mp->prec);
    mpfi_set_q(mp->mval, mp->qval);
    mpq_clear(mp->qval);
    mp->qsize = 0;
}

static void mpqi_canonicalize(mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpqi_arg_check(mp);
    else {
        mpfr_t val;
        mpfr_init2(val, mp->prec);
        mpfi_diam_abs(val, mp->mval);
        if (mpfr_sgn(val) == 0) {
            mpfi_get_left(val, mp->mval);
            mpq_init(mp->qval);
            mpfr_get_q(mp->qval, val);
            mp->qsize = mpq_bytes(mp->qval);
            if (size_ok(mp->qsize)) {
                mpfi_clear(mp->mval);
            } else {
                mpq_clear(mp->qval);
                mp->qsize = 0;
            }
        }
        mpfr_clear(val);
    }
    opcount++;
}

void mpqi_set_d(mpqi_ptr mp, double d) {
    if (mp->qsize == 0) {
        mpfi_clear(mp->mval);
        mpq_init(mp->qval);
    }
    mpq_set_d(mp->qval, d);
    mp->qsize = mpq_bytes(mp->qval);
    mpqi_canonicalize(mp);
}

void mpqi_set_q(mpqi_ptr mp, mpq_srcptr q) {
    if (mp->qsize == 0) {
        mpfi_clear(mp->mval);
        mpq_init(mp->qval);
    }
    mpq_set(mp->qval, q);
    mp->qsize = mpq_bytes(mp->qval);
    mpqi_canonicalize(mp);
}

void mpqi_set_m(mpqi_ptr mp, mpfi_srcptr m) {
    if (mp->qsize > 0) {
        mpq_clear(mp->qval);
        mpfi_init2(mp->mval, mp->prec);
    }
    mpfi_set(mp->mval, m);
    mp->qsize = 0;
    mpqi_canonicalize(mp);
}

void mpqi_set(mpqi_ptr mp, mpqi_ptr m) {
    if (m->qsize > 0)
        mpqi_set_q(mp, m->qval);
    else
        mpqi_set_m(mp, m->mval);
}

double mpqi_get_d(mpqi_ptr mp) {
    if (mp->qsize > 0)
        return mpq_get_d(mp->qval);
    else
        return mpfi_get_d(mp->mval);
}

void mpqi_left(mpfr_ptr left, mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpfr_set_q(left, mp->qval, MPFR_RNDD);
    else
        mpfi_get_left(left, mp->mval);
}

void mpqi_right(mpfr_ptr right, mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpfr_set_q(right, mp->qval, MPFR_RNDU);
    else
        mpfi_get_right(right, mp->mval);
}


void mpqi_mid(mpfr_ptr mid, mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpfr_set_q(mid, mp->qval, MPFR_RNDN);
    else
        mpfi_mid(mid, mp->mval);
}

void mpqi_mid_q(mpq_ptr mid, mpqi_ptr mp) {
    if (mp->qsize > 0)
        mpq_set(mid, mp->qval);
    else {
        mpfr_t mmid;
        mpfr_init2(mmid, mp->prec);
        mpfi_mid(mmid, mp->mval);
        mpfr_get_q(mid, mmid);
        mpfr_clear(mmid);
    }
}

void mpqi_mul_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q) {
    mpqi_arg_check(arg);
    bool clear_m = false;
    bool clear_q = false;
    if (arg->qsize > 0) {
        if (dest->qsize == 0) {
            clear_m = true;
            mpq_init(dest->qval);
        }
        mpq_mul(dest->qval, arg->qval, q);
        dest->qsize = mpq_bytes(dest->qval);
    } else {
        if (dest->qsize > 0) {
            clear_q = true;
            mpfi_init2(dest->mval, dest->prec);
        }
        mpfi_mul_q(dest->mval, arg->mval, q);
        dest->qsize = 0;
    }
    if (clear_m)
        mpfi_clear(dest->mval);
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

void mpqi_mul_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr v) {
    mpqi_arg_check(arg);
    bool clear_q = false;
    if (dest->qsize > 0) {
        clear_q = true;
        mpfi_init2(dest->mval, dest->prec);
    }
    if (arg->qsize > 0)
        mpfi_mul_q(dest->mval, v, arg->qval);
    else
        mpfi_mul(dest->mval, v, arg->mval);
    dest->qsize = 0;
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

void mpqi_mul(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2) {
    mpqi_arg_check(arg2);
    if (arg2->qsize > 0)
        mpqi_mul_q(dest, arg1, arg2->qval);
    else
        mpqi_mul_mpfi(dest, arg1, arg2->mval);
}

void mpqi_add_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q) {
    mpqi_arg_check(arg);
    bool clear_q = false;
    bool clear_m = false;
    if (arg->qsize > 0) {
        if (dest->qsize == 0) {
            clear_m = true;
            mpq_init(dest->qval);
        }
        mpq_add(dest->qval, arg->qval, q);
        dest->qsize = mpq_bytes(dest->qval);
    } else {
        if (dest->qsize > 0) {
            clear_q = true;
            mpfi_init2(dest->mval, dest->prec);
        }
        mpfi_add_q(dest->mval, arg->mval, q);
        dest->qsize = 0;
    }
    if (clear_m)
        mpfi_clear(dest->mval);
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

void mpqi_add_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr v) {
    mpqi_arg_check(arg);
    bool clear_q = false;
    if (dest->qsize > 0) {
        clear_q = true;
        mpfi_init2(dest->mval, dest->prec);
    }
    if (arg->qsize > 0)
        mpfi_add_q(dest->mval, v, arg->qval);
    else
        mpfi_add(dest->mval, v, arg->mval);
    dest->qsize = 0;
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

void mpqi_add(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2) {
    mpqi_arg_check(arg2);
    if (arg2->qsize > 0)
        mpqi_add_q(dest, arg1, arg2->qval);
    else
        mpqi_add_mpfi(dest, arg1, arg2->mval);
}

// written by AI -- unverified
void mpqi_sub_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q) {
    mpqi_arg_check(arg);
    bool clear_m = false;
    bool clear_q = false;
    if (arg->qsize > 0) {
        if (dest->qsize == 0) {
            clear_m = true;
            mpq_init(dest->qval);
        }
        mpq_sub(dest->qval, arg->qval, q);
        dest->qsize = mpq_bytes(dest->qval);
    } else {
        if (dest->qsize > 0) {
            clear_q = true;
            mpfi_init2(dest->mval, dest->prec);
        }
        mpfi_sub_q(dest->mval, arg->mval, q);
        dest->qsize = 0;
    }
    if (clear_m)
        mpfi_clear(dest->mval);
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

// written by AI -- unverified
void mpqi_sub_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr v) {
    mpqi_arg_check(arg);
    bool clear_q = false;
    if (dest->qsize > 0) {
        clear_q = true;
        mpfi_init2(dest->mval, dest->prec);
    }
    if (arg->qsize > 0) {
        mpfi_set_q(dest->mval, arg->qval);
        mpfi_sub(dest->mval, dest->mval, v);
    } else
        mpfi_sub(dest->mval, arg->mval, v);
    dest->qsize = 0;
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

// written by AI -- unverified
void mpqi_sub(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2) {
    mpqi_arg_check(arg2);
    if (arg2->qsize > 0)
        mpqi_sub_q(dest, arg1, arg2->qval);
    else
        mpqi_sub_mpfi(dest, arg1, arg2->mval);
}

// written by AI -- unverified
void mpqi_div_q(mpqi_ptr dest, mpqi_ptr arg, mpq_srcptr q) {
    if (mpq_sgn(q) == 0) {
        mpqi_arg_check(arg);
        bool had_qval = false;
        if (dest->qsize > 0) {
            had_qval = true;
            mpfi_init2(dest->mval, dest->prec);
        }
        mpfr_t left, right;
        mpfr_init2(left, dest->prec);
        mpfr_init2(right, dest->prec);
        mpfr_set_inf(left, -1);
        mpfr_set_inf(right, +1);
        mpfi_interv_fr(dest->mval, left, right);
        mpfr_clear(left);
        mpfr_clear(right);
        dest->qsize = 0;
        if (had_qval)
            mpq_clear(dest->qval);
        mpqi_canonicalize(dest);
        return;
    }

    mpq_t inv;
    mpq_init(inv);
    mpq_inv(inv, q);
    mpqi_mul_q(dest, arg, inv);
    mpq_clear(inv);
}

// written by AI -- unverified
void mpqi_div_mpfi(mpqi_ptr dest, mpqi_ptr arg, mpfi_srcptr v) {
    mpqi_arg_check(arg);
    bool clear_q = false;
    if (dest->qsize > 0) {
        clear_q = true;
        mpfi_init2(dest->mval, dest->prec);
    }
    if (arg->qsize > 0) {
        mpfi_set_q(dest->mval, arg->qval);
        mpfi_div(dest->mval, dest->mval, v);
    } else
        mpfi_div(dest->mval, arg->mval, v);
    dest->qsize = 0;
    if (clear_q)
        mpq_clear(dest->qval);
    mpqi_canonicalize(dest);
}

// written by AI -- unverified
void mpqi_div(mpqi_ptr dest, mpqi_ptr arg1, mpqi_ptr arg2) {
    mpqi_arg_check(arg2);
    if (arg2->qsize > 0)
        mpqi_div_q(dest, arg1, arg2->qval);
    else
        mpqi_div_mpfi(dest, arg1, arg2->mval);
}

// written by AI -- unverified
bool mpqi_is_zero(mpqi_ptr mp) {
    if (mp->qsize > 0)
        return mpq_sgn(mp->qval) == 0;
    else
        return mpfi_is_zero(mp->mval) != 0;
}

// written by AI -- unverified
bool mpqi_has_zero(mpqi_ptr mp) {
    if (mp->qsize > 0)
        return mpq_sgn(mp->qval) == 0;
    else
        return mpfi_has_zero(mp->mval) != 0;
}

double digit_precision_mpqi(mpqi_ptr mp) {
    if (mp->qsize > 0)
        return MAX_DIGIT_PRECISION;
    mpfr_prec_t prec = mpfi_get_prec(mp->mval);
    mpfr_t left;
    mpfr_init2(left, prec);
    mpfr_t right;
    mpfr_init2(right, prec);
    mpfi_get_left(left, mp->mval);
    mpfi_get_right(right, mp->mval);
    if (mpfr_sgn(left) != mpfr_sgn(right))
        return 0.0;
    mpfr_t diam;
    mpfr_init2(diam, prec);
    mpfi_diam_rel(diam, mp->mval);
    if (mpfr_sgn(diam) == 0) {
        mpfr_clear(diam);
        return MAX_DIGIT_PRECISION;
    }
    mpfr_log10(diam, diam, MPFR_RNDN);
    double result = -mpfr_get_d(diam, MPFR_RNDN);
    if (result <= 0)
        result = 0.0;
    if (result > MAX_DIGIT_PRECISION)
        result = MAX_DIGIT_PRECISION;
    mpfr_clear(diam);
    return result;
}
