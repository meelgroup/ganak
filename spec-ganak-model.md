# The ganak-model language

`ganak-model` reads a probabilistic model (`.gm`), translates it into a weighted,
projected CNF and, with `--solve`, computes its queries with Ganak. The web page
compiles the same CNF to a d-DNNF and answers every query on the circuit.

```
# Flight-control hydraulics, per flight. Illustrative numbers.
event eng1  "Engine 1 shuts down"   2e-5
event dual  "Both engines lost"     1e-7
eng1_out = eng1 | dual
green "Green pressure" = !leak_g & (g_src | y_src & !ptu)
lost  "Loss of control" = !(elevator & ailerons & rudder)
query lost given eng1_out
require lost < 1e-9
```

## 1. Lexical structure

- `#` starts a comment that runs to the end of the line.
- A statement ends at a newline or `;`. A newline does not end it inside `( )`
  or `{ }`, or after a token that cannot end an expression (a binary operator,
  `,`, `=`, `!`, `given`, `then`, `else`, `if`, `of`).
- Names: `[A-Za-z_][A-Za-z0-9_]*`. A choice value is `choice.value`.
- Labels: `"..."`, with `\"` and `\\` escapes. Labels are for display only.
- Probabilities: `0.25`, `.25`, `2e-5`, `1/3`, `25%`. They are read exactly
  (`2e-5` is 1/50000) and must lie in [0, 1].
- Reserved words: `event choice rule constraint query require given if then else
  true false not and or xor flip atleast atmost exactly of`.

## 2. Statements

| Statement | Meaning |
|---|---|
| `event NAME [LABEL] PROB` | an independent Boolean random event, true with probability PROB |
| `choice NAME [LABEL] { V1 P1, V2 P2, ... }` | an independent random variable with values V1, V2, ...; the Pi must sum to exactly 1; `NAME.Vi` is true iff the value is Vi |
| `NAME [LABEL] = EXPR` | a definition |
| `rule BODY -> NAME [PROB]` | BODY causes NAME with probability PROB (default 1) |
| `constraint EXPR` | a hard constraint: every probability is conditioned on it |
| `query EXPR [given E1, E2, ...]` | report P(EXPR \| E1 ∧ E2 ∧ ...) |
| `require EXPR [given E1, ...] OP PROB` | as `query`, and check it: OP is `<`, `<=`, `>`, `>=` |

Names are global and can be used before they are defined. A name is declared
once, except that a name with rules may also be declared an `event` (see 3.2).

## 3. Expressions

From loosest to tightest binding:

| Operator | Meaning | Associativity |
|---|---|---|
| `if C then A else B` | (C ∧ A) ∨ (¬C ∧ B); the `else` branch extends as far as possible | |
| `<->` | equivalence | left |
| `->` | implication | right |
| `\|`, `or` | disjunction | left |
| `^`, `xor` | exclusive or | left |
| `&`, `and` | conjunction | left |
| `!`, `~`, `not` | negation | prefix |

Atoms: a name, `true`, `false`, `( EXPR )`, `flip(PROB)` (a fresh independent
event, one per textual occurrence), and `atleast K [of] (E1, ..., En)`,
`atmost K [of] (...)`, `exactly K [of] (...)`.

### 3.1 The probability space

The random variables are the events, the choices, every `flip`, every rule with
a probability below 1, and the spontaneous part of a rule head that is also an
event. They are mutually independent. Definitions and rules must not form a
cycle, so every assignment to the random variables determines every name. Each
assignment has its product probability P₀. With C the conjunction of all
constraints:

    P(Q | E) = P₀(Q ∧ E ∧ C) / P₀(E ∧ C)

The answer is undefined when the denominator is 0, and is reported as such.
`ganak-model` also reports P₀(C), the probability that the constraints hold.

### 3.2 Rules

The rules `rule B1 -> h p1`, ..., `rule Bn -> h pn` define

    h = (B1 & flip(p1)) | ... | (Bn & flip(pn))      [ | h.self ]

with `flip(1)` read as `true`. If `event h p` is also declared, its own
randomness is the event `h.self` with probability p, which joins the
disjunction ("leak"). This is ProbLog's semantics for acyclic propositional
programs, and noisy-OR in Bayesian-network terms. A name cannot have both rules
and a `=` definition.

## 4. The CNF

The output is a weighted projected CNF in the model counting competition format
(`c t pwmc`):

- One variable per random event, per choice value, per name, and per query and
  `given` expression that is not already a literal. Helper variables of
  subexpressions are shared by structural hashing.
- Every definition is encoded in both directions (full Tseitin), so every
  variable is a function of the random variables.
- An event with probability p: `c p weight v p 0` and `c p weight -v 1-p 0`.
  A choice value with probability p: `c p weight v p 0`, `c p weight -v 1 0`,
  plus clauses saying exactly one value holds.
- `c p show` lists the random variables, the names, and the query variables, so
  `bdnnf --assume` can condition on any of them. Helper variables are not shown.
- `c <var> <name>` names each shown variable (anonymous ones are named
  `flip@LINE:COL`, `rule@LINE`).
- `c gm query <i> <lit> [given <lits>]` and
  `c gm require <i> <lit> [given <lits>] <op> <bound>` record the statements.
- Constraints are clauses (unit clauses on their literals).

The weighted projected count is P₀(C). P(Q | E) is the count under
assumptions Q and E divided by the count under E.

## 5. Command line

```
ganak-model model.gm                 # the CNF on stdout
ganak-model model.gm -o model.cnf --json model.json --dot model.dot
ganak-model --solve model.gm         # every query and requirement, exact
```

Errors are printed as `FILE:LINE:COL: error: ...`, and the exit code is 1.

`--json` writes the model for the web page: the nodes (name, label, kind,
probability, variable, the names they depend on, source line), the statements,
and the errors.
