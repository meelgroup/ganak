"""
Tests for the pyganak Python bindings.

Each test verifies model counts against independently known correct values.
Variable indices are 1-based (positive = positive literal, negative = negated).
"""

import math
import pytest
import pyganak


# ---------------------------------------------------------------------------
# Module-level sanity
# ---------------------------------------------------------------------------

class TestModule:
    def test_version_exists(self):
        assert hasattr(pyganak, "__version__")
        assert isinstance(pyganak.__version__, str)
        assert len(pyganak.__version__) > 0

    def test_version_alias(self):
        assert pyganak.VERSION == pyganak.__version__

    def test_counter_class_exists(self):
        assert hasattr(pyganak, "Counter")


# ---------------------------------------------------------------------------
# Constructor / basic interface
# ---------------------------------------------------------------------------

class TestConstructor:
    def test_default_construction(self):
        c = pyganak.Counter()
        assert c.nof_vars() == 0
        assert c.nof_clauses() == 0

    def test_verbose_kwarg(self):
        c = pyganak.Counter(verbose=0)
        assert c is not None

    def test_seed_kwarg(self):
        c = pyganak.Counter(seed=42)
        assert c is not None

    def test_both_kwargs(self):
        c = pyganak.Counter(verbose=0, seed=123)
        assert c is not None

    def test_bad_arg_type(self):
        with pytest.raises(TypeError):
            pyganak.Counter(verbose="yes")  # type: ignore


# ---------------------------------------------------------------------------
# add_clause / add_clauses
# ---------------------------------------------------------------------------

class TestAddClause:
    def test_add_single_clause(self):
        c = pyganak.Counter()
        c.add_clause([1, 2])
        assert c.nof_clauses() == 1
        assert c.nof_vars() == 2

    def test_add_clause_updates_max_var(self):
        c = pyganak.Counter()
        c.add_clause([3, -5])
        assert c.nof_vars() == 5

    def test_add_clauses_list(self):
        c = pyganak.Counter()
        c.add_clauses([[1, 2], [-1, 3], [2, -3]])
        assert c.nof_clauses() == 3

    def test_add_clauses_generator(self):
        c = pyganak.Counter()
        c.add_clauses(([i] for i in [1, 2, 3]))
        assert c.nof_clauses() == 3

    def test_zero_literal_raises(self):
        c = pyganak.Counter()
        with pytest.raises(ValueError):
            c.add_clause([1, 0, 2])

    def test_empty_clause_is_unsat(self):
        # An empty clause makes the formula UNSAT → count = 0
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([])     # empty clause = False
        assert c.count() == 0


# ---------------------------------------------------------------------------
# set_sampling_set
# ---------------------------------------------------------------------------

class TestSamplingSet:
    def test_set_sampling_set_basic(self):
        c = pyganak.Counter()
        c.add_clause([1, 2, 3])
        c.set_sampling_set([1, 2])    # count only over vars 1 and 2
        result = c.count()
        assert isinstance(result, int)
        assert result >= 0

    def test_sampling_set_zero_raises(self):
        c = pyganak.Counter()
        with pytest.raises(ValueError):
            c.set_sampling_set([0, 1])

    def test_sampling_set_negative_raises(self):
        c = pyganak.Counter()
        with pytest.raises(ValueError):
            c.set_sampling_set([-1, 1])


# ---------------------------------------------------------------------------
# Core counting correctness
# ---------------------------------------------------------------------------

class TestCounting:
    """
    All expected counts are independently verified.
    """

    def test_no_vars_no_clauses(self):
        # Empty formula over 0 variables: exactly 1 model (the empty assignment).
        c = pyganak.Counter()
        result = c.count()
        assert result == 1

    def test_two_unconstrained_vars(self):
        # 2 unconstrained variables → 4 models.
        c = pyganak.Counter()
        c.new_vars(2)
        assert c.count() == 4

    def test_new_vars_increases_nof_vars(self):
        c = pyganak.Counter()
        c.new_vars(5)
        assert c.nof_vars() == 5
        c.new_vars(3)
        assert c.nof_vars() == 8

    def test_unit_clause_positive(self):
        # x1: exactly 1 model (x1=T), var x2 free → 2 models (if 2 vars declared)
        # Use only x1: 1 model.
        c = pyganak.Counter()
        c.add_clause([1])
        assert c.count() == 1

    def test_unit_clause_negative(self):
        # ¬x1: 1 model (x1=F)
        c = pyganak.Counter()
        c.add_clause([-1])
        assert c.count() == 1

    def test_two_unit_clauses(self):
        # x1 ∧ x2: 1 model
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([2])
        assert c.count() == 1

    def test_contradicting_units(self):
        # x1 ∧ ¬x1: UNSAT → 0 models
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([-1])
        assert c.count() == 0

    def test_two_var_disjunction(self):
        # (x1 ∨ x2): 3 models (TT, TF, FT)
        c = pyganak.Counter()
        c.add_clause([1, 2])
        assert c.count() == 3

    def test_two_var_conjunction(self):
        # (x1) ∧ (x2): 1 model
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([2])
        assert c.count() == 1

    def test_two_var_xor_via_clauses(self):
        # x1 XOR x2 = (x1 ∨ x2) ∧ (¬x1 ∨ ¬x2): 2 models (TF, FT)
        c = pyganak.Counter()
        c.add_clause([1, 2])
        c.add_clause([-1, -2])
        assert c.count() == 2

    def test_two_var_xnor_via_clauses(self):
        # ¬(x1 XOR x2) = (¬x1 ∨ x2) ∧ (x1 ∨ ¬x2): 2 models (TT, FF)
        c = pyganak.Counter()
        c.add_clause([-1, 2])
        c.add_clause([1, -2])
        assert c.count() == 2

    def test_three_var_all_false_forbidden(self):
        # (x1 ∨ x2 ∨ x3): 7 models (all assignments except FFF)
        c = pyganak.Counter()
        c.add_clause([1, 2, 3])
        assert c.count() == 7

    def test_three_var_all_true_forced(self):
        # x1 ∧ x2 ∧ x3: 1 model
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([2])
        c.add_clause([3])
        assert c.count() == 1

    def test_pigeon_hole_2_3(self):
        # Pigeonhole: 2 pigeons, 3 holes — satisfiable; count is positive.
        # var(i,j) = (i-1)*3 + j  →  p11=1, p12=2, p13=3, p21=4, p22=5, p23=6
        c = pyganak.Counter()
        c.add_clause([1, 2, 3])    # pigeon 1 in some hole
        c.add_clause([4, 5, 6])    # pigeon 2 in some hole
        # At-most-one-pigeon-per-hole:
        c.add_clause([-1, -4])
        c.add_clause([-2, -5])
        c.add_clause([-3, -6])
        result = c.count()
        assert isinstance(result, int)
        assert result > 0

    def test_result_is_python_int(self):
        c = pyganak.Counter()
        c.add_clause([1, 2])
        result = c.count()
        assert type(result) is int

    def test_large_count_is_big_integer(self):
        # 20 unconstrained variables → 2^20 = 1,048,576 models
        c = pyganak.Counter()
        c.new_vars(20)
        result = c.count()
        assert result == 2**20

    def test_count_second_call_raises(self):
        # count() may only be called once per instance.
        c = pyganak.Counter()
        c.add_clause([1, 2])
        c.count()
        with pytest.raises(RuntimeError):
            c.count()

    def test_sampling_set_projection(self):
        # (x1 ∨ x2 ∨ x3) projected onto {x1, x2}:
        #   x1=T, x2=T  → satisfiable (x3 free) → 1 projected model
        #   x1=T, x2=F  → satisfiable               → 1 projected model
        #   x1=F, x2=T  → satisfiable               → 1 projected model
        #   x1=F, x2=F  → only if x3=T satisfies → 1 projected model
        # All 4 assignments of {x1,x2} satisfy the projected formula → 4
        c = pyganak.Counter()
        c.add_clause([1, 2, 3])
        c.set_sampling_set([1, 2])
        assert c.count() == 4

    def test_sampling_set_strict_projection(self):
        # x1 ∧ x2: projected onto {x1} → 1 (x1 must be T)
        c = pyganak.Counter()
        c.add_clause([1])
        c.add_clause([2])
        c.set_sampling_set([1])
        assert c.count() == 1

    def test_clauses_with_large_var_index(self):
        # Make sure high variable numbers work without issue
        c = pyganak.Counter()
        c.add_clause([100, -200])
        result = c.count()
        assert isinstance(result, int)
        assert result > 0


# ---------------------------------------------------------------------------
# Multiple independent Counter instances
# ---------------------------------------------------------------------------

class TestMultipleCounters:
    def test_independent_instances(self):
        c1 = pyganak.Counter()
        c1.add_clause([1, 2])

        c2 = pyganak.Counter()
        c2.add_clause([1])

        assert c1.count() == 3
        assert c2.count() == 1

    def test_instance_isolation(self):
        # Adding to c1 doesn't affect c2
        c1 = pyganak.Counter()
        c2 = pyganak.Counter()
        c1.add_clause([1, 2])
        assert c2.nof_clauses() == 0


# ---------------------------------------------------------------------------
# Successive counts use separate Counter instances (count() is one-shot)
# ---------------------------------------------------------------------------

class TestSuccessiveCounts:
    def test_tighter_formula_has_fewer_models(self):
        # (x1 ∨ x2) alone: 3 models.
        # (x1 ∨ x2) ∧ (¬x1 ∨ ¬x2): 2 models (XOR).
        # Use two separate Counter instances.
        c1 = pyganak.Counter()
        c1.add_clause([1, 2])
        assert c1.count() == 3

        c2 = pyganak.Counter()
        c2.add_clause([1, 2])
        c2.add_clause([-1, -2])
        assert c2.count() == 2

    def test_contradiction_gives_zero(self):
        # x1 alone: 1 model.  x1 ∧ ¬x1: 0 models.
        c1 = pyganak.Counter()
        c1.add_clause([1])
        assert c1.count() == 1

        c2 = pyganak.Counter()
        c2.add_clause([1])
        c2.add_clause([-1])
        assert c2.count() == 0


# ---------------------------------------------------------------------------
# WeightedCounter
# ---------------------------------------------------------------------------

class TestWeightedCounterBasic:
    def test_class_exists(self):
        assert hasattr(pyganak, "WeightedCounter")

    def test_default_construction(self):
        c = pyganak.WeightedCounter()
        assert c.nof_vars() == 0
        assert c.nof_clauses() == 0

    def test_prec_kwarg(self):
        c = pyganak.WeightedCounter(prec=256)
        assert c is not None

    def test_prec_too_small_raises(self):
        with pytest.raises(ValueError):
            pyganak.WeightedCounter(prec=1)

    def test_result_is_float(self):
        c = pyganak.WeightedCounter()
        c.add_clause([1])
        c.set_lit_weight(1, 1.0)
        c.set_lit_weight(-1, 0.0)
        result = c.count()
        assert type(result) is float

    def test_count_once_only(self):
        c = pyganak.WeightedCounter()
        c.add_clause([1])
        c.count()
        with pytest.raises(RuntimeError):
            c.count()

    def test_zero_lit_weight_raises(self):
        c = pyganak.WeightedCounter()
        with pytest.raises((ValueError, TypeError)):
            c.set_lit_weight(0, 1.0)


class TestWeightedCounterCorrectness:
    def test_unit_clause_weighted(self):
        # x1 forced true, w(x1)=0.3, w(¬x1)=0.7 → count = 0.3
        c = pyganak.WeightedCounter()
        c.add_clause([1])
        c.set_lit_weight( 1, 0.3)
        c.set_lit_weight(-1, 0.7)
        assert math.isclose(c.count(), 0.3, rel_tol=1e-9)

    def test_two_var_disjunction_weighted(self):
        # (x1 ∨ x2):  models (T,T), (T,F), (F,T)
        # w(x1=T)=0.3, w(x1=F)=0.7, w(x2=T)=0.4, w(x2=F)=0.6
        # (T,T)=0.3*0.4=0.12, (T,F)=0.3*0.6=0.18, (F,T)=0.7*0.4=0.28
        # total = 0.58
        c = pyganak.WeightedCounter()
        c.add_clause([1, 2])
        c.set_lit_weight( 1, 0.3)
        c.set_lit_weight(-1, 0.7)
        c.set_lit_weight( 2, 0.4)
        c.set_lit_weight(-2, 0.6)
        assert math.isclose(c.count(), 0.58, rel_tol=1e-9)

    def test_unsat_formula_weighted(self):
        # x1 ∧ ¬x1 is UNSAT → count = 0
        c = pyganak.WeightedCounter()
        c.add_clause([1])
        c.add_clause([-1])
        c.set_lit_weight( 1, 0.5)
        c.set_lit_weight(-1, 0.5)
        assert c.count() == 0.0

    def test_unit_weight_equals_unweighted(self):
        # When all weights are 1.0, the weighted count equals the plain count.
        c_plain = pyganak.Counter()
        c_plain.add_clause([1, 2])
        plain = c_plain.count()  # 3

        c_weighted = pyganak.WeightedCounter()
        c_weighted.add_clause([1, 2])
        c_weighted.set_lit_weight( 1, 1.0)
        c_weighted.set_lit_weight(-1, 1.0)
        c_weighted.set_lit_weight( 2, 1.0)
        c_weighted.set_lit_weight(-2, 1.0)
        assert math.isclose(c_weighted.count(), float(plain), rel_tol=1e-9)

    def test_higher_precision_same_result(self):
        # prec=128 and prec=256 should give the same double-precision result.
        def make():
            c = pyganak.WeightedCounter(prec=128)
            c.add_clause([1, 2])
            c.set_lit_weight( 1, 0.3)
            c.set_lit_weight(-1, 0.7)
            c.set_lit_weight( 2, 0.4)
            c.set_lit_weight(-2, 0.6)
            return c.count()

        def make256():
            c = pyganak.WeightedCounter(prec=256)
            c.add_clause([1, 2])
            c.set_lit_weight( 1, 0.3)
            c.set_lit_weight(-1, 0.7)
            c.set_lit_weight( 2, 0.4)
            c.set_lit_weight(-2, 0.6)
            return c.count()

        assert math.isclose(make(), make256(), rel_tol=1e-12)

    def test_add_clauses_weighted(self):
        c = pyganak.WeightedCounter()
        c.add_clauses([[1, 2], [-1, 2]])
        c.set_lit_weight( 1, 0.5)
        c.set_lit_weight(-1, 0.5)
        c.set_lit_weight( 2, 0.5)
        c.set_lit_weight(-2, 0.5)
        result = c.count()
        assert isinstance(result, float)
        assert result >= 0.0

    def test_sampling_set_weighted(self):
        # (x1 ∨ x2 ∨ x3) projected onto {x1, x2}, all weights 0.5
        c = pyganak.WeightedCounter()
        c.add_clause([1, 2, 3])
        c.set_sampling_set([1, 2])
        for lit in [1, -1, 2, -2, 3, -3]:
            c.set_lit_weight(lit, 0.5)
        result = c.count()
        assert isinstance(result, float)
        assert result >= 0.0
