// Unit tests for the interval Fields (FMpfi, FMpqi): is_zero/is_one must be exact.

#include "../src/fmpqi.hpp"

#include <cstdio>

static int failures = 0;
static int checks = 0;

#define CHECK(cond) do {                                                    \
    ++checks;                                                               \
    if (!(cond)) {                                                          \
        fprintf(stderr, "  FAIL [line %d]: %s\n", __LINE__, #cond);        \
        ++failures;                                                         \
    }                                                                       \
} while (0)

static const mpfr_prec_t prec = 64;

// [lo, lo + 2^-exp]
static void set_interval(mpfi_t m, long lo, long exp) {
    mpfr_t l, r;
    mpfr_init2(l, prec);
    mpfr_init2(r, prec);
    mpfr_set_si(l, lo, MPFR_RNDN);
    mpfr_set_si_2exp(r, 1, -exp, MPFR_RNDN);
    mpfr_add(r, r, l, MPFR_RNDU);
    mpfi_interv_fr(m, l, r);
    mpfr_clear(l);
    mpfr_clear(r);
}

static std::unique_ptr<FMpfi> mpfi_interval(long lo, long exp) {
    auto f = std::make_unique<FMpfi>(prec);
    set_interval(f->val, lo, exp);
    return f;
}

static std::unique_ptr<FMpqi> mpqi_interval(long lo, long exp) {
    mpfi_t m;
    mpfi_init2(m, prec);
    set_interval(m, lo, exp);
    auto f = std::make_unique<FMpqi>(prec);
    mpqi_set_m(&f->val, m);
    mpfi_clear(m);
    return f;
}

static void test_fmpfi() {
    printf("Testing: FMpfi is_zero/is_one\n");
    CHECK(FMpfi(0L, prec).is_zero());
    CHECK(!FMpfi(0L, prec).is_one());
    CHECK(FMpfi(1L, prec).is_one());
    CHECK(!FMpfi(1L, prec).is_zero());
    CHECK(!FMpfi(2L, prec).is_one());

    CHECK(!mpfi_interval(1, 60)->is_one());
    CHECK(!mpfi_interval(0, 60)->is_one());
    CHECK(!mpfi_interval(0, 60)->is_zero());
    CHECK(!mpfi_interval(-1, 0)->is_zero());  // [-1, 0]
    FMpfi around_one(prec);
    mpfi_interv_d(around_one.val, 0.5, 1.5);
    CHECK(!around_one.is_one());
    FMpfi around_zero(prec);
    mpfi_interv_d(around_zero.val, -0.5, 0.5);
    CHECK(!around_zero.is_zero());
}

static void test_fmpqi() {
    printf("Testing: FMpqi is_zero/is_one\n");
    CHECK(FMpqi(0L, prec).is_zero());
    CHECK(!FMpqi(0L, prec).is_one());
    CHECK(FMpqi(1L, prec).is_one());
    CHECK(!FMpqi(1L, prec).is_zero());
    CHECK(!FMpqi(2L, prec).is_one());

    auto near_one = mpqi_interval(1, 60);
    CHECK(near_one->val.qsize == 0);
    CHECK(!near_one->is_one());
    auto near_zero = mpqi_interval(0, 60);
    CHECK(near_zero->val.qsize == 0);
    CHECK(!near_zero->is_zero());
    CHECK(!near_zero->is_one());
}

int main() {
    printf("=== interval Field unit tests ===\n\n");
    test_fmpfi();
    test_fmpqi();
    printf("\n=== Results: %d checks, %d failures ===\n", checks, failures);
    return failures > 0 ? 1 : 0;
}
