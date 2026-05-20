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
 * Pure-C equivalent of example.cpp:
 *   - 10 variables, clause (x1 v x2)
 *   - weighted: positive literal of x4 has weight 10, negative has weight 1
 *   - sampling/independent set = {x1, x2, x4}
 *   - run Arjun preprocessing, then run Ganak, print the count.
 *
 * Expected output: "count is: 33"
 */

#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

#include <arjun/arjun_c.h>
#include <ganak/ganak_c.h>

static void die(const char* what) {
    fprintf(stderr, "error: %s: %s\n", what, arjun_last_error());
    fprintf(stderr, "       (ganak: %s)\n", ganak_last_error());
    exit(1);
}

int main(void) {
    /* ---------------- Build a SimplifiedCNF over rationals -------------- */
    arjun_fgen_t* fg = arjun_fgen_mpq_new();
    if (!fg) die("arjun_fgen_mpq_new");

    arjun_simpcnf_t* cnf = arjun_simpcnf_new(fg);
    if (!cnf) die("arjun_simpcnf_new");

    arjun_simpcnf_new_vars(cnf, 10);

    /* Clause: (x1 v x2)  -- DIMACS-style signed ints. */
    int32_t cl[] = {1, 2};
    if (arjun_simpcnf_add_clause(cnf, cl, 2) != 0) die("add_clause");

    /* Weighted: w(x4=T) = 10, w(x4=F) = 1. */
    arjun_simpcnf_set_weighted(cnf, 1);

    arjun_field_t* w_pos = arjun_field_from_int(fg, 10);
    arjun_field_t* w_neg = arjun_field_from_int(fg, 1);
    if (!w_pos || !w_neg) die("field_from_int");
    if (arjun_simpcnf_set_lit_weight(cnf,  4, w_pos) != 0) die("set weight +4");
    if (arjun_simpcnf_set_lit_weight(cnf, -4, w_neg) != 0) die("set weight -4");
    arjun_field_free(w_pos);
    arjun_field_free(w_neg);

    /* Sampling vars (0-indexed): x1->0, x2->1, x4->3. */
    uint32_t sampl[] = {0, 1, 3};
    if (arjun_simpcnf_set_sampl_vars(cnf, sampl, 3) != 0) die("set sampl");

    /* ---------------- Run Arjun preprocessing --------------------------- */
    arjun_arjun_t* arj = arjun_new();
    if (!arj) die("arjun_new");
    arjun_set_verb(arj, 0);
    if (arjun_standalone_minimize_indep(arj, cnf, /*all_indep=*/0) != 0)
        die("minimize_indep");
    if (arjun_standalone_elim_to_file(arj, cnf) != 0) die("elim_to_file");

    /* ---------------- Run Ganak ----------------------------------------- */
    ganak_counter_t* counter = ganak_counter_new(/*verb=*/0, /*seed=*/0, fg);
    if (!counter) die("ganak_counter_new");
    if (ganak_counter_setup_from_simpcnf(counter, cnf) != 0) die("setup");
    arjun_field_t* sub_count = ganak_counter_count(counter,
                                                   /*bits_threads=*/0,
                                                   /*num_threads=*/1,
                                                   /*debug_threads=*/0);
    if (!sub_count) die("count");

    /* ---------------- Combine with the Arjun multiplier ----------------- */
    arjun_field_t* total = arjun_simpcnf_get_multiplier_weight(cnf);
    if (!total) die("get_multiplier_weight");
    if (!arjun_field_is_zero(total)) {
        if (arjun_field_mul_assign(total, sub_count) != 0) die("mul");
    }

    /* ---------------- Print the count ----------------------------------- */
    mpq_t out;
    mpq_init(out);
    if (arjun_field_get_mpq(total, out) != 0) die("get_mpq");
    gmp_printf("count is: %Qd\n", out);
    mpq_clear(out);

    /* ---------------- Cleanup ------------------------------------------- */
    arjun_field_free(total);
    arjun_field_free(sub_count);
    ganak_counter_free(counter);
    arjun_free(arj);
    arjun_simpcnf_free(cnf);
    arjun_fgen_free(fg);
    return 0;
}
