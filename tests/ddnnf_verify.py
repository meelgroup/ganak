#!/usr/bin/env python3
"""Verify a d4 .nnf circuit from `ganak --compile`.

parse() -> (nodes, arcs, root); count() = structural count (OR=sum, AND=product,
t=1, f=0, arc lits ignored); models() = exact model set (small inputs) using arc
lits. Counting is structural: free vars are explicit OR(v,-v), so no smoothing.
"""
import sys


def parse(path):
    nodes = {}                  # id -> 'o'|'a'|'t'|'f'
    arcs = {}                   # id -> list of (child_id, [lits])
    root = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            toks = line.split()
            if toks[0] == 'c':
                continue
            if toks[0] in ('o', 'a', 't', 'f'):
                nid = int(toks[1])
                nodes[nid] = toks[0]
                arcs.setdefault(nid, [])
                if root is None:
                    root = nid      # first declared node is the root
            else:
                parent = int(toks[0])
                child = int(toks[1])
                lits = [int(x) for x in toks[2:] if x != '0']
                arcs.setdefault(parent, []).append((child, lits))
    return nodes, arcs, root


def count(nodes, arcs, root):
    memo = {}

    def rec(nid):
        if nid in memo:
            return memo[nid]
        t = nodes[nid]
        if t == 't':
            r = 1
        elif t == 'f':
            r = 0
        elif t == 'o':
            r = sum(rec(c) for c, _ in arcs[nid])
        else:  # 'a'
            r = 1
            for c, _ in arcs[nid]:
                r *= rec(c)
        memo[nid] = r
        return r

    return rec(root)


def models(nodes, arcs, root, nvars):
    """Return set of full models (each a frozenset of signed dimacs literals).
    Only call on small circuits/instances."""
    memo = {}

    def merge(a, b):
        # a, b: dict var->bool ; return merged or None on conflict
        out = dict(a)
        for v, val in b.items():
            if v in out and out[v] != val:
                return None
            out[v] = val
        return out

    def apply_lits(m, lits):
        out = dict(m)
        for l in lits:
            v = abs(l)
            val = l > 0
            if v in out and out[v] != val:
                return None
            out[v] = val
        return out

    def rec(nid):
        if nid in memo:
            return memo[nid]
        t = nodes[nid]
        if t == 't':
            res = [dict()]
        elif t == 'f':
            res = []
        elif t == 'o':
            res = []
            for c, lits in arcs[nid]:
                for m in rec(c):
                    mm = apply_lits(m, lits)
                    if mm is not None:
                        res.append(mm)
        else:  # 'a'
            res = [dict()]
            for c, lits in arcs[nid]:
                child_models = []
                for m in rec(c):
                    mm = apply_lits(m, lits)
                    if mm is not None:
                        child_models.append(mm)
                new = []
                for base in res:
                    for cm in child_models:
                        merged = merge(base, cm)
                        if merged is not None:
                            new.append(merged)
                res = new
        memo[nid] = res
        return res

    partial = rec(root)
    full = set()
    for m in partial:
        free = [v for v in range(1, nvars + 1) if v not in m]
        for mask in range(1 << len(free)):
            assign = dict(m)
            for i, v in enumerate(free):
                assign[v] = bool((mask >> i) & 1)
            full.add(frozenset(v if assign[v] else -v for v in range(1, nvars + 1)))
    return full


def evaluate(nodes, arcs, root, assign):
    """Evaluate the circuit as a Boolean function on a complete assignment
    (var->bool, 1-indexed); an arc fires only if its literals hold. The semantics
    functional synthesis needs; correct regardless of decomposability."""
    memo = {}

    def lit_true(l):
        return (l > 0) == assign[abs(l)]

    def rec(nid):
        if nid in memo:
            return memo[nid]
        t = nodes[nid]
        if t == 't':
            r = True
        elif t == 'f':
            r = False
        elif t == 'o':
            r = any(all(lit_true(l) for l in lits) and rec(c) for c, lits in arcs[nid])
        else:  # 'a'
            r = all(all(lit_true(l) for l in lits) and rec(c) for c, lits in arcs[nid])
        memo[nid] = r
        return r

    return rec(root)


def function_models(nodes, arcs, root, nvars):
    """Models of the circuit as a function (brute force over all assignments)."""
    out = set()
    for bits in __import__('itertools').product((False, True), repeat=nvars):
        assign = {v + 1: bits[v] for v in range(nvars)}
        if evaluate(nodes, arcs, root, assign):
            out.add(frozenset(v if assign[v] else -v for v in range(1, nvars + 1)))
    return out


def synthesize(nodes, arcs, root, xassign):
    """Functional synthesis: given input assignment `xassign` (var->bool over X),
    follow the circuit to a satisfying path and return a full assignment, reading
    the Y vars off arc literals. None if no path is consistent with xassign."""

    def lit_ok(l):
        v = abs(l)
        if v in xassign:        # input variable: must match X
            return (l > 0) == xassign[v]
        return True             # output variable: free to take the arc's value

    def add(a, lits):
        for l in lits:
            v, b = abs(l), l > 0
            if v in a and a[v] != b:
                return False
            a[v] = b
        return True

    def rec(nid):
        t = nodes[nid]
        if t == 't':
            return {}
        if t == 'f':
            return None
        if t == 'o':
            for c, lits in arcs[nid]:
                if not all(lit_ok(l) for l in lits):
                    continue
                sub = rec(c)
                if sub is None:
                    continue
                a = dict(sub)
                if add(a, lits):
                    return a
            return None
        # AND: all children must agree
        a = {}
        for c, lits in arcs[nid]:
            if not all(lit_ok(l) for l in lits):
                return None
            sub = rec(c)
            if sub is None or not add(a, lits):
                return None
            for v, b in sub.items():
                if v in a and a[v] != b:
                    return None
                a[v] = b
        return a

    return rec(root)


def subtree_vars(nodes, arcs, root):
    """node -> (pos_vars, neg_vars) appearing on arcs in its reachable subgraph."""
    memo = {}

    def rec(nid):
        if nid in memo:
            return memo[nid]
        memo[nid] = (set(), set())   # guard against cycles (there are none, but safe)
        pos, neg = set(), set()
        for c, lits in arcs.get(nid, []):
            for l in lits:
                (pos if l > 0 else neg).add(abs(l))
            cp, cn = rec(c)
            pos |= cp
            neg |= cn
        memo[nid] = (pos, neg)
        return memo[nid]

    rec(root)
    return memo


def reachable_nodes(nodes, arcs, root):
    """Set of node ids reachable from the root."""
    seen = set()
    stack = [root]
    while stack:
        nid = stack.pop()
        if nid in seen:
            continue
        seen.add(nid)
        for c, _ in arcs.get(nid, []):
            stack.append(c)
    return seen


def unreachable_nodes(nodes, arcs, root):
    """Declared-but-unreachable ("dead") node ids. Empty for a cleaned circuit;
    raw `ganak --compile` output may have some. Asserts ddnnf-cleanup drops them."""
    return set(nodes) - reachable_nodes(nodes, arcs, root)


def check_decomposable(nodes, arcs, root):
    """Strict d-DNNF: every AND node's children have pairwise-disjoint variable sets."""
    sv = subtree_vars(nodes, arcs, root)
    for nid, t in nodes.items():
        if t != 'a':
            continue
        if nid not in sv:           # unreachable ("dead") node: not part of the circuit
            continue
        seen = set()
        for c, _ in arcs.get(nid, []):
            cp, cn = sv[c]
            cvars = cp | cn
            if seen & cvars:
                return False, f"AND {nid} children share var(s) {sorted(seen & cvars)}"
            seen |= cvars
    return True, "ok"


def check_weak_decomposable(nodes, arcs, root):
    """Weak (Akshay) d-DNNF: no variable appears positively in one child of an AND
    node and negatively in another child."""
    sv = subtree_vars(nodes, arcs, root)
    for nid, t in nodes.items():
        if t != 'a':
            continue
        if nid not in sv:           # unreachable ("dead") node: not part of the circuit
            continue
        kids = [sv[c] for c, _ in arcs.get(nid, [])]
        for i in range(len(kids)):
            for j in range(len(kids)):
                if i == j:
                    continue
                bad = kids[i][0] & kids[j][1]   # pos in i, neg in j
                if bad:
                    return False, f"AND {nid}: var(s) {sorted(bad)} both polarities across children"
    return True, "ok"


def shared_and_vars(nodes, arcs, root):
    """Return the set of variables that appear in more than one child subtree of
    some AND node (i.e. the variables on which decomposability is relaxed)."""
    sv = subtree_vars(nodes, arcs, root)
    shared = set()
    for nid, t in nodes.items():
        if t != 'a':
            continue
        if nid not in sv:           # unreachable ("dead") node: not part of the circuit
            continue
        seen = set()
        for c, _ in arcs.get(nid, []):
            cp, cn = sv[c]
            cvars = cp | cn
            shared |= (seen & cvars)
            seen |= cvars
    return shared


def read_cnf(path):
    """Parse a small DIMACS CNF -> (nvars, clauses). Comment/sampling lines
    ignored; only for brute-forcing tiny unprojected fixtures."""
    nvars = 0
    clauses = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or line[0] == 'c':
                continue
            if line[0] == 'p':
                nvars = int(line.split()[2])
                continue
            lits = [int(x) for x in line.split() if x not in ('', '0')]
            if lits:
                clauses.append(lits)
    return nvars, clauses


def brute_models(nvars, clauses):
    """Exact model set of a CNF by brute force (frozenset of signed lits per
    model). Only for tiny nvars."""
    import itertools
    out = set()
    for bits in itertools.product((False, True), repeat=nvars):
        assign = {v + 1: bits[v] for v in range(nvars)}
        if all(any((l > 0) == assign[abs(l)] for l in cl) for cl in clauses):
            out.add(frozenset(v if assign[v] else -v for v in range(1, nvars + 1)))
    return out


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(
        description="Verify a d4 .nnf circuit from `ganak --compile`. With no "
                    "check flag, just prints the structural model count.")
    ap.add_argument('nnf')
    ap.add_argument('--expect-count', type=int, default=None,
                    help="assert the circuit's structural model count equals N")
    ap.add_argument('--cnf', default=None,
                    help="brute-force this small UNPROJECTED CNF and assert the "
                         "circuit's count AND full model set match it exactly "
                         "(also validates arc literals)")
    ap.add_argument('--check-decomposable', action='store_true',
                    help="assert strict decomposability (AND children var-disjoint)")
    ap.add_argument('--check-weak-decomposable', action='store_true',
                    help="assert weak decomposability (no var both polarities across AND children)")
    ap.add_argument('--strict', action='store_true',
                    help="assert canonical cleaned form: no dead nodes, root id 1, ids 1..N contiguous")
    ap.add_argument('-q', '--quiet', action='store_true')
    args = ap.parse_args()

    nodes, arcs, root = parse(args.nnf)
    c = count(nodes, arcs, root)

    any_check = (args.expect_count is not None or args.cnf or args.check_decomposable
                 or args.check_weak_decomposable or args.strict)
    if not any_check:
        print(c)               # backward-compatible: bare invocation prints the count
        sys.exit(0)

    def fail(msg):
        print("ddnnf-verify: FAIL: " + msg)
        sys.exit(1)

    if args.expect_count is not None and c != args.expect_count:
        fail("count %d != expected %d" % (c, args.expect_count))

    if args.cnf:
        nvars, clauses = read_cnf(args.cnf)
        true_models = brute_models(nvars, clauses)
        if c != len(true_models):
            fail("count %d != #models(CNF) %d" % (c, len(true_models)))
        circ = models(nodes, arcs, root, nvars)
        if circ != true_models:
            fail("model set mismatch: +%d extra, -%d missing"
                 % (len(circ - true_models), len(true_models - circ)))

    if args.check_decomposable:
        ok, msg = check_decomposable(nodes, arcs, root)
        if not ok:
            fail("not decomposable: " + msg)

    if args.check_weak_decomposable:
        ok, msg = check_weak_decomposable(nodes, arcs, root)
        if not ok:
            fail("not weak-decomposable: " + msg)

    if args.strict:
        dead = unreachable_nodes(nodes, arcs, root)
        if dead:
            fail("%d dead node(s) present: %s" % (len(dead), sorted(dead)[:8]))
        if root != 1:
            fail("root id %d != 1" % root)
        if set(nodes) != set(range(1, len(nodes) + 1)):
            fail("node ids not contiguous 1..N")

    if not args.quiet:
        print("ddnnf-verify: OK (count=%d, nodes=%d)" % (c, len(nodes)))
    sys.exit(0)
