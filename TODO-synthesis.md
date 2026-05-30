# TODO: `--synthesis` — component overlap over inputs (X)

Status as of 2026-05-30. Context for this file: `--synthesis 1` does Boolean
functional synthesis via wDNNF share-and-branch. See the synthesis notes in
`CLAUDE.md` and the README `--synthesis` blurb for the current mechanism and the
two reference papers.

## Background: what the papers permit

`--synthesis` implements Boolean functional synthesis in the sense of:

- **Akshay, Arora, Chakraborty, Krishna, Raghunathan, Shah**, *Knowledge
  Compilation for Boolean Functional Synthesis*, FMCAD 2019 —
  <https://www.cse.iitb.ac.in/~supratik/publications/papers/FMCAD19.pdf>
  (introduces **SynNNF**; recaps **wDNNF** from Akshay et al. 2018).
- **Illner, Kučera**, *A Compiler for Weak Decomposable Negation Normal Form*,
  AAAI 2024 — <https://ojs.aaai.org/index.php/AAAI/article/download/28926/29761>.

Two definitions matter:

- **wDNNF** (Akshay et al. 2018 / AAAI24 Def. 1): an AND node `α = α₁∧…∧αₖ` is
  *weak decomposable* if every variable shared by two children is **pure** — only
  `x` or only `¬x` appears in the subcircuit under `α`
  (`lits(αᵣ) ∩ {¬ℓ | ℓ ∈ lits(αₛ)} = ∅`). **No input/output distinction** — purely
  a polarity condition.
- **SynNNF** (FMCAD19, the actual poly-time synthesis condition, weaker than
  wDNNF): the `∧ᵢ`-unrealizability constraint is imposed **only on the synthesized
  (output) variables**. The inputs carry **no decomposability constraint at all**.

> ⚠️ Convention clash: FMCAD19 writes `F(X, Y)` with `X` = **outputs**, `Y` =
> **inputs** — the *opposite* of Ganak's code, where `X` = inputs/sampling vars
> (`< indep_support_end`) and `Y` = outputs (`>= indep_support_end`).

## The question

Current code shares **only output vars** (`>= indep_support_end`), keeping inputs
disjoint (`CompAnalyzer::compute_shareable_vars`, the `v < indep_end` /
`v >= indep_end` guards at `src/comp_analyzer.cpp:256,263,270`). Can we also overlap
components over **inputs**?

## Verdict

**Theory: yes.** Neither wDNNF nor SynNNF requires inputs to stay disjoint
(wDNNF allows sharing any pure var; SynNNF constrains only the outputs). Ganak's
"only outputs may be shared" rule is **strictly more conservative than the theory
demands.**

**But it is not a one-line guard removal**, because of an asymmetry the abstract
definitions hide — it's about *how Ganak realizes "sharing"*, not the normal form:

- Ganak's "share" = mark the var **non-bridging** + **pin it to its pure polarity**
  (`Counter::synth_forced_lit`). Operationally this is **"fix the variable to one
  value."**
- Fixing is sound only for variables we are **free to choose** — the
  existentially-quantified **outputs** (the SAT oracle supplies the `Y` witness;
  pinning makes sibling SAT leaves agree).
- **Inputs are universally given by the synthesis query.** The circuit must
  faithfully represent **both** polarities of every input. Pinning/fixing a shared
  input deletes the `x=false` branch → the circuit becomes unfaithful as a function
  of the inputs. **So the output mechanism cannot be reused for inputs.**

The genuinely-wDNNF way to share an input is **branch-not-fix**: keep it a real
branch in each sibling (same polarity, never fixed). That *is* faithful for
synthesis — `φ₁|ₓ ∧ φ₂|ₓ = F|ₓ` for either value of `x`, since conjunction
distributes — and only breaks the count (already meaningless in this mode). This
faithfulness argument does not even need purity, matching SynNNF's "no constraint
on inputs."

**Conclusion:** the restriction is over-conservative w.r.t. the papers, but
*correct* given Ganak's fix-by-pinning realization of sharing. Lifting it requires
a **second, branch-not-fix sharing path** distinct from the output path.

### Caveats that bound the payoff

1. Ganak branches the **independent support (inputs) first**, so at most decompose
   points the inputs are already assigned — few unknown inputs remain to share. The
   succinctness upside is likely modest. Measure before investing.

## Plan

1. **Fix the existing soundness bug first (prerequisite).** ✅ DONE. The fix:
   - `compute_shareable_vars` refuses to share an output var whose original CNF
     has BOTH polarities (`orig_polarity[v] == 3`) — no globally-consistent pin
     exists for such a var.
   - `synth_forced_lit` falls back to `orig_polarity[v]` (state-independent)
     when the residual is FREE, rather than defaulting to `+v`.
   - `ddnnf.hpp::check_wdnnf_at_and` (SLOW_DEBUG only) asserts the wDNNF
     invariant at every `mk_and` to catch regressions immediately.
   `tests/ddnnf_synth.py --synthesis` is at 0 failures across thousands of
   random projected instances; ctest `ddnnf_synth_share_and_branch` exercises it.
2. **Add a gated flag** (e.g. `--synthesis 2`, or `--synthshareinputs 1`) so
   input-sharing is opt-in and A/B-testable against current output-only behaviour.
3. **Branch-not-fix path in the analyzer.** In `compute_shareable_vars`, allow
   unknown input vars (`v < indep_end`) to be marked shareable, but tag them as a
   **distinct kind** ("share-by-branch") vs the output "share-by-fix" kind. Keep the
   demotion fixpoint as-is (it is var-agnostic; ensures every active clause keeps a
   bridging owner so nothing is dropped).
4. **Suppress pinning for shared inputs.** `synth_forced_lit` must return
   `lit_Undef` for share-by-branch (input) vars — they get branched normally in each
   sibling, never forced. Pinning stays output-only.
5. **Confirm the trace mechanics** allow a shared input to be branched independently
   in each sibling: it is already re-marked unvisited in `make_comp_from_archetype`,
   but verify backtracking re-frees it cleanly between siblings and the component
   cache is not keyed in a way that conflates the two sub-searches.
6. **Validate empirically** (decisive step, per repo policy — fuzz after every
   change):
   - `tests/ddnnf_synth.py --synthesis --num 200` — round-trip soundness; compare
     the circuit size ratio with/without input-sharing to judge whether the payoff
     justifies the complexity.
   - `tests/ddnnf_fuzz.py 200` — faithful d-DNNF not regressed.
   - `../count_fuzzer/fuzz.py --only 200 --exact` — no count regressions in the
     shared analyzer code paths.

## Relevant code

- `src/comp_analyzer.cpp:242` — `compute_shareable_vars` (the `indep_end` guards).
- `src/comp_analyzer.cpp:336` — `record_comp` (`is_shareable` → non-bridging).
- `src/comp_analyzer.hpp:200,210` — share-mode hooks; re-mark shared vars unvisited.
- `src/counter.hpp:204`, `src/counter.cpp` — `synth_forced_lit` (pure-polarity pin).
- `src/comp_management.hpp:75` — `comp_has_shareable`.
- `src/counter_config.hpp:43` — `synthesis` config field.
</content>
</invoke>
