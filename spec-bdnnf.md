# BDNNF: Binary d-DNNF format, version 1

BDNNF (`.bdnnf`) stores a d-DNNF (or Decision-DNNF) circuit. It is designed so
that:

- a compiler can write it in **one pass, append-only, never seeking back**,
  without knowing the number of nodes, edges or the root in advance (it can be
  written to a pipe);
- a consumer can use it **without loading it into memory**: bottom-up queries
  (counting, weighted counting, counting under assumptions) are one forward
  sequential scan, top-down queries are one backward sequential scan, and any
  node can be located through a small block index;
- it is compact: node ids are implicit, child references are relative,
  literals are gap-encoded varints, and blocks are optionally zstd-compressed.

## 1. Logical model

Nodes are numbered `0, 1, 2, ...`. Nodes `0` and `1` are the constants FALSE
and TRUE; they are never stored. Stored nodes are numbered `2, 3, ...` in the
order they appear in the file. **Every child of a node has a smaller id than
the node itself**, so the file is in topological order, children first.

Variables are `1..nvars`. A literal is a non-zero integer: `v` is variable `v`
true, `-v` is variable `v` false.

There are two kinds of stored nodes.

**AND**: `AND(children; units; free)`, where `children` are node ids, `units`
are literals and `free` are variables:

    f = (AND of f(c) over children) AND (AND of units) AND (AND of (v OR -v) over free)

A free variable is logically a tautology, but it is part of the node's scope:
it is how a compiler says "this variable is unconstrained here", which matters
for counting (factor 2 per free variable).

**OR**: `OR(arcs)`, where each arc is a pair `(lits, child)`:

    f = OR over arcs of ((AND of lits) AND f(child))

In a Decision-DNNF an OR has two arcs: one whose literals contain `x` and one
whose literals contain `-x`, for the decision variable `x`; the other literals
on each arc are those implied by the decision. The format does not require
this shape.

The **root** is given in the footer (section 4). It can be any node, including
`0` or `1` (e.g. an unsatisfiable formula has root `0` and no stored nodes).
Nodes not reachable from the root ("dead nodes") are allowed; a writer can
avoid declaring them by flagging the file `CLEAN` (section 4.1).

### 1.1 Defined variables

A file may also carry **definitions** (section 3.1): a function, given as an
AIG, for each of some variables. A defined variable is not mentioned by any
node. The models of the file are the models `m` of the circuit, each extended
by setting every defined variable `v` to its function `f_v(m)`; a function may
use other defined variables, acyclically. Counting queries are answered by the
circuit alone, as long as no defined variable is counted over: a compiler uses
definitions for the variables outside the projection set (Skolem functions of
the variables it eliminated), so a projected count, conditioning on the
projection set, and sampling all stay plain circuit queries, and a sampled
model is completed by the definitions.

### 1.2 Structural count

    count(FALSE) = 0,  count(TRUE) = 1
    count(AND)   = product of count(c) over children  *  2^|free|
    count(OR)    = sum of count(child) over arcs

Literals are ignored. This is the model count over `nvars` variables when the
file is flagged `DECOMPOSABLE`, `DETERMINISTIC` and `SMOOTH` (section 4.1).
Weighted counting is the same with `count(AND)` additionally multiplied by
the weight of each unit, `(w(v) + w(-v))` for each free variable, and the
weights of the arc literals multiplied into each OR arc.

## 2. File layout

All fixed-width integers are **little-endian**. `varint` is unsigned LEB128
(7 bits per byte, least significant group first, high bit set on all but the
last byte), at most 10 bytes for 64 bits.

    Header
    Frame*         node frames, then exactly one footer frame
    Trailer        16 bytes

### 2.1 Header

| Bytes | Field |
|---|---|
| 8 | magic `89 42 44 4E 4E 46 0D 0A` (`\x89BDNNF\r\n`) |
| 4 | `u32 version` = 1 |
| 4 | `u32 hlen`: length of the header payload that follows |
| hlen | header payload |

Header payload:

    varint nvars
    varint nkv
    nkv * (varint klen, klen bytes key, varint vlen, vlen bytes value)

Keys and values are UTF-8 strings. Defined keys: `producer` (e.g.
`ganak 1.2.3 <sha>`), and `projection`: the projection variables, decimal,
space-separated (may be empty; absent: all of `1..nvars`). Projected queries (count, sample,
conditioning) are over these variables. Readers ignore unknown keys. Ganak
also writes `input` (the CNF file), `command` (its command line), and the
version (SHA1) of each component: `arjun`, `sbva`, `cryptominisat`, `cadical`,
`cadiback`, `approxmc`, `treedecomp`. Unknown `version` values
must be rejected.

The non-ASCII first byte and the `\r\n` let tools detect text-mode
corruption, as in PNG.

### 2.2 Frames

Every frame has a 20-byte head, a payload, and a 4-byte tail:

| Bytes | Field |
|---|---|
| 1 | `u8 kind`: 1 = node block, 2 = footer, 3 = definitions |
| 1 | `u8 codec`: 0 = raw, 1 = zstd (one zstd frame) |
| 2 | `u16 reserved` = 0 |
| 4 | `u32 nnodes`: number of node records in this frame (0 for the footer) |
| 4 | `u32 stored_len`: length of the stored payload |
| 4 | `u32 raw_len`: length of the payload after decompression |
| 4 | `u32 crc`: CRC-32C (Castagnoli) of the *raw* (decompressed) payload |
| stored_len | payload |
| 4 | `u32 stored_len` again |

The repeated `stored_len` lets a reader walk the frames backwards from the
footer to the first frame. The total frame size is `stored_len + 24`.

A node record is never split across frames. A writer typically closes a frame
when its raw payload reaches 64–256 KiB. The first node record of the first
node frame is node `2`; the first node of each later frame follows the last
node of the previous one.

## 3. Node records

A node block's raw payload is a sequence of `nnodes` records. Each record
starts with a tag byte:

| Tag bits | Meaning |
|---|---|
| 0–1 | kind: 0 = AND, 1 = OR (2, 3 reserved) |
| 2 | AND only: has units |
| 3 | AND only: has free variables |
| 4–7 | reserved, must be 0 |

**AND record**

    tag
    varint nkids, nkids * ref
    if tag bit 2: litlist     (units)
    if tag bit 3: varlist     (free variables)

**OR record**

    tag
    varint narcs, narcs * (litlist, ref)

**ref** is a child reference relative to the current node id `self`:

    varint r:  r = 0  -> FALSE (node 0)
               r = 1  -> TRUE  (node 1)
               r >= 2 -> node  self - r + 1, which must be in [2, self-1]

so `r = 2` is the previous node. Compilers that work depth-first refer mostly
to recent nodes, so most refs fit in one or two bytes.

**litlist**: `varint n`, then the `n` literals sorted by their code
`c(l) = 2*|l| + (l < 0 ? 1 : 0)`, strictly increasing: the first code is
written as is, each later one as `c_i - c_{i-1} - 1`.

**varlist**: `varint n`, then the `n` variables strictly increasing: the
first as is, each later one as `v_i - v_{i-1} - 1`.

### 3.1 Definitions frame

At most one, right before the footer. Its raw payload:

    varint nnodes                AND nodes of the AIG
    nnodes * (aref, aref)        the two inputs of each AND node
    varint ndefs
    ndefs * (varint var, aref)   var := aref

An **aref** is a varint `r`: `k = r >> 1` is TRUE for `k = 0`, variable `k` for
`1 <= k <= nvars`, and the `j`-th AND node for `k = nvars + 1 + j`; the
reference is negated if `r & 1`. An AND node refers only to variables and to
earlier AND nodes. Each variable is defined at most once.

## 4. Footer

The footer frame's raw payload:

    varint nnodes        stored nodes (so the largest id is nnodes + 1)
    varint nedges        AND children + OR arcs
    varint nlits         literals in units and arc litlists, plus free variables
    varint root          root node id (0, 1, or 2..nnodes+1)
    varint props         property flags, section 4.1
    varint nframes       node frames
    nframes * (varint frame_size, varint frame_nnodes)
    varint defs_size     full size of the definitions frame, 0 if there is none
    varint nkv
    nkv * (varint klen, key, varint vlen, value)

`frame_size` is the full size of each node frame (`stored_len + 24`); the
first node frame starts right after the header. Summing sizes gives the file
offset of each frame, and summing `frame_nnodes` gives the first node id in
each frame, so a reader can find the frame of any node by binary search.

Defined footer keys: `count` (the model count the writer computed, decimal;
lets a consumer check its own count). Readers ignore unknown keys.

### 4.1 Property flags

| Bit | Name | Promise |
|---|---|---|
| 0 | `DECOMPOSABLE` | For every AND reachable from the root, its children, its units and its free variables mention pairwise disjoint variable sets; and for every OR arc, its literals are on distinct variables that its child does not mention. |
| 1 | `DETERMINISTIC` | For every OR reachable from the root, the arcs are pairwise logically inconsistent. |
| 2 | `SMOOTH` | Every reachable OR's arcs mention the same variables, and the root mentions all of `1..nvars`. Together with bits 0–1 this makes the structural count the model count. |
| 3 | `CLEAN` | No dead nodes, and the root is the last stored node. |

A writer sets only the bits it guarantees. Other bits are reserved and are 0.

## 5. Trailer

| Bytes | Field |
|---|---|
| 8 | `u64` file offset of the footer frame |
| 8 | magic `42 44 4E 4E 46 45 4E 44` (`BDNNFEND`) |

## 6. Reading

- **Streaming (pipe)**: read the header, then frames until the footer frame.
  Everything needed is in order; the root arrives last, with the footer.
- **Random access**: read the trailer, then the footer, then seek to frames
  using the frame index. Decoded frames can be kept in an LRU cache.
- **Backward scan**: from the footer frame, read the 4 bytes before each
  frame's head to find the previous frame's `stored_len`.
- **No trailer** (the writer was killed): scan frames from the header; the
  CRCs show where a partially written frame begins. Such a file has no root
  and is incomplete, but the frames before the cut can be verified.

Memory for a bottom-up evaluation is one value per node. A post-processing
tool can lower this by recording, for each node, its last use (the largest id
that references it), so the evaluator can drop values once they are dead.

## 7. Text form

`bdnnf dump` prints a file as text, one node per line, with ids absolute and
`F`/`T` for the constants:

    bdnnf 1 nvars 3
    2 and free 2
    3 or 1 -> 2 ; -1 2 -> T
    4 and 3 free 3
    root 4
    props DECOMPOSABLE DETERMINISTIC SMOOTH CLEAN
    meta count 6

This is `(x1 OR x2)` over 3 variables. Node 3 decides on x1: with x1 true, x2
is free (node 2, count 2); with x1 false, x2 is implied (count 1). Node 4
conjoins that with the free variable 3: (2 + 1) * 2 = 6 models.

AND lines are `<id> and <kids...> [units <lits...>] [free <vars...>]`. OR
lines are `<id> or <arc> ; <arc> ...`, where an arc is `<lits...> -> <child>`.
Header and footer key/value pairs are printed as `meta <key> <value>` lines,
before the nodes and after `props` respectively. Definitions are printed as
`defnode <k> <aref> <aref>` (the AND node referred to as `k`) and
`def <var> <aref>` lines, before `root`.
