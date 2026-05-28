# Component cache vs. share-and-branch wDNNF (`--weak 3`)

Why normal component caching is unsound for the synthesis share-and-branch
mode, what invariant it breaks, and what a sound fix can and cannot recover.

## The cache invariant

The component cache key is

    (unassigned member vars) + (active / unsatisfied clause-ids)

and its soundness rests on a single property:

> **The key uniquely determines the residual sub-formula over those vars.**

Two occurrences of the same key are then literally the same sub-problem, so the
cached count (and, in compile mode, the cached d-DNNF node) is interchangeable.

In a variable-disjoint decomposition this property holds because of two facts
that work together:

1. **Connected closure.** Every active clause touching a member has *all* of its
   unassigned literals inside the member set — component analysis bridges
   through every variable until closure. So there are no dangling unassigned
   literals outside the key.
2. **Assigned ⟹ false.** Any *assigned* literal in an active clause must be
   false (otherwise the clause would be satisfied and excluded). So the assigned
   context surrounding the component is *implied* by the key, not free.

Together: same key ⟹ same residual problem ⟹ interchangeable result.

## How `--weak 3` breaks it

Share-and-branch deliberately lets an input variable belong to several sibling
components without bridging. In `CompAnalyzer::record_comp`:

    if (is_shareable(v)) { claimed_share[v] = 1; continue; }

a shared input is recorded as a *member* but we do **not** bridge through it.
That kills invariant #1 for shared vars: the clauses that constrain a shared
member can live in a *sibling* component and are therefore **not in this
component's key**.

The break shows up via **unit propagation**, not via the component's own clause
set. Concretely (a traced 11-var, k=9 case, Y = {10,11}):

- Clause `{2,6,7}`. In a context with inputs `2=F, 7=F`, BCP propagates `6=T`.
  That implied literal `6=T` is baked into the emitted node ("input 7=F ⟹ 6=T")
  and the node is **cached**.
- The variable that forced `6=T` is the shared var **2** — whose clause was
  never bridged into this component, so `2` is *not* in this component's key.
- The node is later reused under `2=T, 7=F`: now `{2,6,7}` is satisfied and `6`
  is free, but the cached node still asserts `6=T`. The corresponding satisfiable
  X (with `6=F`) loses its witness → **under-coverage**.

So the cached node carries an implied literal produced by a clause that is not
in its key. The key no longer determines the residual sub-formula ⟹ reuse is
unsound. Empirically: full caching gives 10–14 no-witness synthesis failures per
200 random projected CNFs.

## The uncomfortable part: the size win *is* the bug

With full caching the circuit is actually **smaller** than `--weak 0`
(~0.96–0.99×). But that win is *produced by* the unsound reuse: reusing one
`{6,7}` node across both the `2=F` and `2=T` leaves is simultaneously the
compaction and the bug. Any keying that makes reuse sound must stop reusing
across differing shared-var contexts — which removes exactly the reuse that
shrank the circuit. You cannot keep *that specific* win by refining the key.

## Fix options

1. **Selective caching (shipped default).** Do not cache any component that
   contains a shared member (`CompManager::comp_has_shareable`, guard in
   `Counter::backtrack`). Sound. Shared components are recomputed each time, so
   the circuit ends up ~1.05–1.10× *larger* than `--weak 0`. Simple and correct.

2. **Context-extended key (the principled refinement).** Fold the *values of the
   assigned shared vars adjacent to the component* into the fingerprint — i.e.
   the external context that can produce implied literals. Then a node built
   under `2=F` is keyed differently from a query under `2=T` (sound), but it
   still hits across all the other inputs whose value doesn't touch this
   component, so most of the reuse survives (≈2 variants of the `{6,7}` node
   instead of 1, not one-per-leaf). Sound and plausibly still compact. Cost:
   during analysis, collect assigned shared vars adjacent to a member and fold
   them with polarity into the `DifferencePackedComp` fingerprint — medium
   effort, touches the packed-comp format. This is the only sound option that
   can plausibly beat `--weak 0`.

3. **Never collapse a shared var ("true wDNNF").** Force each shared member to be
   branched *inside* every component that owns it, so the node is a genuine
   function of it and is sound to reuse for any value. Already measured to run
   ~1.05× *bigger* (duplicating the shared subtree + losing its caching), so it
   does not pay.

## Recommendation

Keep selective caching (option 1) as the shipped sound default. If the
compaction matters, prototype the adjacent-assigned-shared-var key extension
(option 2) and measure size + soundness on `tests/ddnnf_synth.py`.
