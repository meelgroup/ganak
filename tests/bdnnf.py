#!/usr/bin/env python3
"""Independent Python reader, checker and counter for BDNNF (spec-bdnnf.md).

Usage: bdnnf.py check FILE | count FILE [--cnf CNF] [--structural] | dump FILE
"""
import struct
import sys
from fractions import Fraction

MAGIC = b"\x89BDNNF\r\n"
END_MAGIC = b"BDNNFEND"
PROPS = ["DECOMPOSABLE", "DETERMINISTIC", "SMOOTH", "CLEAN"]
P_DECOMPOSABLE, P_DETERMINISTIC, P_SMOOTH, P_CLEAN = 1, 2, 4, 8


class FormatError(Exception):
    pass


def _crc_table():
    t = []
    for i in range(256):
        c = i
        for _ in range(8):
            c = (c >> 1) ^ (0x82F63B78 if c & 1 else 0)
        t.append(c)
    return t


_CRC = _crc_table()


def crc32c(data):
    c = 0xFFFFFFFF
    t = _CRC
    for b in data:
        c = (c >> 8) ^ t[(c ^ b) & 0xFF]
    return c ^ 0xFFFFFFFF


def _zstd_decompress(data):
    try:
        from compression import zstd
        return zstd.decompress(data)
    except ImportError:
        pass
    try:
        import zstandard
        return zstandard.ZstdDecompressor().decompress(data)
    except ImportError:
        raise FormatError("zstd frame, but no zstd module (need Python 3.14 or 'zstandard')")


class Buf:
    def __init__(self, data):
        self.d = data
        self.p = 0

    def done(self):
        return self.p == len(self.d)

    def byte(self):
        if self.p >= len(self.d):
            raise FormatError("truncated record")
        b = self.d[self.p]
        self.p += 1
        return b

    def varint(self):
        v = 0
        sh = 0
        while True:
            b = self.byte()
            v |= (b & 0x7F) << sh
            if not b & 0x80:
                break
            sh += 7
            if sh > 63:
                raise FormatError("varint too long")
        if v >= 1 << 64:
            raise FormatError("varint overflows 64 bits")
        return v

    def string(self):
        n = self.varint()
        if self.p + n > len(self.d):
            raise FormatError("truncated string")
        s = self.d[self.p:self.p + n].decode("utf-8")
        self.p += n
        return s

    def kv(self):
        return [(self.string(), self.string()) for _ in range(self.varint())]


class Circuit:
    """nodes[id] for id >= 2: ('and', kids, units, free) or ('or', [(lits, child)])."""

    def __init__(self):
        self.nvars = 0
        self.header_kv = []
        self.footer_kv = []
        self.nodes = [None, None]
        self.root = 0
        self.props = 0
        self.nedges = 0
        self.nlits = 0
        self.frames = []        # (offset, size, first_id, nnodes, codec)

    def meta(self, key):
        for k, v in self.footer_kv:
            if k == key:
                return v
        return None


def _codes(b):
    n = b.varint()
    out = []
    prev = None
    for i in range(n):
        g = b.varint()
        cur = g if i == 0 else prev + g + 1
        out.append(cur)
        prev = cur
    return out


def _lits(b, nvars, nid):
    out = []
    for c in _codes(b):
        v = c >> 1
        if v < 1 or v > nvars:
            raise FormatError("literal out of range at node %d" % nid)
        out.append(-v if c & 1 else v)
    return out


def _ref(b, nid):
    r = b.varint()
    if r <= 1:
        return r
    child = nid - r + 1
    if child < 2:
        raise FormatError("child reference out of range at node %d" % nid)
    return child


def _decode_node(b, nid, nvars):
    tag = b.byte()
    kind = tag & 3
    if kind == 0:
        if tag & ~0xF:
            raise FormatError("reserved tag bits set at node %d" % nid)
        kids = [_ref(b, nid) for _ in range(b.varint())]
        units = _lits(b, nvars, nid) if tag & 4 else []
        free = []
        if tag & 8:
            free = _codes(b)
            for v in free:
                if v < 1 or v > nvars:
                    raise FormatError("free variable out of range at node %d" % nid)
        return ('and', kids, units, free)
    if kind == 1:
        if tag & ~3:
            raise FormatError("reserved tag bits set at node %d" % nid)
        arcs = []
        for _ in range(b.varint()):
            lits = _lits(b, nvars, nid)
            arcs.append((lits, _ref(b, nid)))
        return ('or', arcs)
    raise FormatError("unknown node kind %d at node %d" % (kind, nid))


def read(path):
    """Parses and validates everything the spec makes checkable from the bytes."""
    with open(path, "rb") as f:
        data = f.read()
    c = Circuit()
    if len(data) < 16 or data[:8] != MAGIC:
        raise FormatError("not a BDNNF file (bad magic)")
    ver, hlen = struct.unpack_from("<II", data, 8)
    if ver != 1:
        raise FormatError("unsupported version %d" % ver)
    if 16 + hlen > len(data):
        raise FormatError("truncated header")
    hb = Buf(data[16:16 + hlen])
    c.nvars = hb.varint()
    c.header_kv = hb.kv()
    if not hb.done():
        raise FormatError("trailing bytes in header")
    pos = 16 + hlen
    data_start = pos
    nid = 2
    footer = None
    while footer is None:
        if pos + 20 > len(data):
            raise FormatError("unexpected end of file (no footer: incomplete file?)")
        kind, codec, res, nn, slen, rlen, crc = struct.unpack_from("<BBHIIII", data, pos)
        if kind not in (1, 2):
            raise FormatError("bad frame kind %d at offset %d" % (kind, pos))
        if res != 0:
            raise FormatError("reserved frame bits set")
        if pos + 24 + slen > len(data):
            raise FormatError("truncated frame at offset %d" % pos)
        stored = data[pos + 20:pos + 20 + slen]
        (tail,) = struct.unpack_from("<I", data, pos + 20 + slen)
        if tail != slen:
            raise FormatError("frame tail mismatch at offset %d" % pos)
        if codec == 0:
            if slen != rlen:
                raise FormatError("raw frame length mismatch")
            raw = stored
        elif codec == 1:
            raw = _zstd_decompress(stored)
            if len(raw) != rlen:
                raise FormatError("zstd frame decompressed to the wrong size")
        else:
            raise FormatError("unknown codec %d" % codec)
        if crc32c(raw) != crc:
            raise FormatError("frame CRC mismatch at offset %d" % pos)
        if kind == 2:
            if nn != 0:
                raise FormatError("footer frame with nnodes != 0")
            footer = (raw, pos)
        else:
            if nn == 0:
                raise FormatError("empty node frame")
            c.frames.append((pos, slen + 24, nid, nn, codec))
            b = Buf(raw)
            for _ in range(nn):
                c.nodes.append(_decode_node(b, nid, c.nvars))
                nid += 1
            if not b.done():
                raise FormatError("trailing bytes in node frame")
        pos += 24 + slen
    raw, foff = footer
    b = Buf(raw)
    nnodes = b.varint()
    c.nedges = b.varint()
    c.nlits = b.varint()
    c.root = b.varint()
    c.props = b.varint()
    if c.props & ~0xF:
        raise FormatError("unknown property bits")
    nframes = b.varint()
    index = [(b.varint(), b.varint()) for _ in range(nframes)]
    c.footer_kv = b.kv()
    if not b.done():
        raise FormatError("trailing bytes in footer")
    if nnodes != nid - 2:
        raise FormatError("footer nnodes %d != %d stored nodes" % (nnodes, nid - 2))
    if [(s, n) for (_, s, _, n, _) in c.frames] != index:
        raise FormatError("footer frame index does not match the frames")
    if c.root >= nnodes + 2:
        raise FormatError("root out of range")
    if pos + 16 != len(data):
        raise FormatError("trailer missing or garbage after it")
    toff, = struct.unpack_from("<Q", data, pos)
    if data[pos + 8:pos + 16] != END_MAGIC or toff != foff:
        raise FormatError("bad trailer")
    if data_start + sum(s for (s, _) in index) != foff:
        raise FormatError("frame sizes do not add up to the footer offset")
    edges = lits = 0
    for n in c.nodes[2:]:
        if n[0] == 'and':
            edges += len(n[1])
            lits += len(n[2]) + len(n[3])
        else:
            edges += len(n[1])
            lits += sum(len(l) for l, _ in n[1])
    if (edges, lits) != (c.nedges, c.nlits):
        raise FormatError("footer nedges/nlits %d/%d != counted %d/%d"
                          % (c.nedges, c.nlits, edges, lits))
    return c


def children(n):
    return n[1] if n[0] == 'and' else [ch for _, ch in n[1]]


def reachable(c):
    seen = set()
    stack = [c.root] if c.root >= 2 else []
    while stack:
        x = stack.pop()
        if x in seen:
            continue
        seen.add(x)
        stack.extend(ch for ch in children(c.nodes[x]) if ch >= 2)
    return seen


def node_vars(c):
    """vars[id]: bitmask (int, bit v = var v) of the vars below each node."""
    vs = [0, 0]
    for n in c.nodes[2:]:
        m = 0
        if n[0] == 'and':
            for k in n[1]:
                m |= vs[k]
            for l in n[2]:
                m |= 1 << abs(l)
            for v in n[3]:
                m |= 1 << v
        else:
            for lits, ch in n[1]:
                m |= vs[ch]
                for l in lits:
                    m |= 1 << abs(l)
        vs.append(m)
    return vs


def low_var(m):
    return (m & -m).bit_length() - 1


def check_props(c):
    """Verifies the property flags the file claims. Returns a list of errors."""
    errs = []
    reach = reachable(c)
    vs = node_vars(c)
    if c.props & P_CLEAN:
        if len(c.nodes) > 2 and c.root != len(c.nodes) - 1:
            errs.append("CLEAN but the root is not the last node")
        if len(reach) != len(c.nodes) - 2:
            errs.append("CLEAN but %d dead nodes" % (len(c.nodes) - 2 - len(reach)))
    for nid in sorted(reach):
        n = c.nodes[nid]
        if n[0] == 'and' and c.props & P_DECOMPOSABLE:
            parts = [vs[k] for k in n[1]] + [1 << abs(l) for l in n[2]] + [1 << v for v in n[3]]
            seen = 0
            for p in parts:
                if seen & p:
                    errs.append("DECOMPOSABLE but AND %d shares var %d" % (nid, low_var(seen & p)))
                    break
                seen |= p
        if n[0] == 'or' and c.props & P_DECOMPOSABLE:
            for lits, ch in n[1]:
                lm = 0
                bad = False
                for x in lits:
                    bad |= bool(lm >> abs(x) & 1)
                    lm |= 1 << abs(x)
                if bad or lm & vs[ch]:
                    errs.append("DECOMPOSABLE but OR %d arc %s overlaps its child %d"
                                % (nid, lits, ch))
        if n[0] == 'or':
            live = [(l, ch) for l, ch in n[1] if ch != 0]
            if c.props & P_SMOOTH:
                sets = set()
                for l, ch in live:
                    m = vs[ch]
                    for x in l:
                        m |= 1 << abs(x)
                    sets.add(m)
                if len(sets) > 1:
                    errs.append("SMOOTH but OR %d arcs mention different vars" % nid)
            if c.props & P_DETERMINISTIC:
                for i in range(len(n[1])):
                    for j in range(i + 1, len(n[1])):
                        a, b = set(n[1][i][0]), {-x for x in n[1][j][0]}
                        if not a & b:
                            errs.append("DETERMINISTIC but OR %d arcs %d, %d share no "
                                        "complementary literal" % (nid, i, j))
    if c.props & P_SMOOTH and c.root != 0:
        missing = ((1 << (c.nvars + 1)) - 2) & ~vs[c.root]
        if missing:
            errs.append("SMOOTH but the root misses var %d" % low_var(missing))
    return errs


def structural_count(c):
    val = [0, 1]
    for n in c.nodes[2:]:
        if n[0] == 'and':
            x = 1 << len(n[3])
            for k in n[1]:
                x *= val[k]
        else:
            x = sum(val[ch] for _, ch in n[1])
        val.append(x)
    return val[c.root]


def read_cnf_weights(path):
    """Returns (nvars, weights {lit: Fraction} or None, projection set or None),
    with Arjun's semantics: the first weight given to a var sets that literal to
    w and the other one to 1-w; a later line overrides its own literal."""
    nvars = 0
    weights = None
    proj = None
    with open(path) as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if t[0] == 'p' and len(t) >= 3:
                nvars = int(t[2])
            elif t[:3] == ['c', 'p', 'weight']:
                weights = {} if weights is None else weights
                l = int(t[3])
                w = Fraction(t[4])
                if l not in weights and -l not in weights:
                    weights[-l] = 1 - w
                weights[l] = w
            elif t[:3] == ['c', 'p', 'show']:
                proj = set() if proj is None else proj
                proj.update(int(x) for x in t[3:] if x != '0')
            elif t[:2] == ['c', 'ind']:
                proj = set() if proj is None else proj
                proj.update(int(x) for x in t[2:] if x != '0')
    return nvars, weights, proj


def weighted_count(c, weights=None, proj=None):
    """Exact (weighted) model count over `proj` (default all vars), smoothing on
    the fly from per-node var sets, so it does not rely on the SMOOTH flag. Needs
    a decomposable, deterministic circuit. Literals of vars outside `proj` weigh 1."""
    P = set(range(1, c.nvars + 1)) if proj is None else set(proj)
    one = Fraction(1) if weights else 1

    def w(l):
        if abs(l) not in P:
            return one
        if weights is None:
            return one
        return weights.get(l, one)

    def both(v):
        return w(v) + w(-v)

    # A logically FALSE subtree mentions no vars: counting its vars into an OR's
    # scope would make that scope overlap AND siblings further up. Scopes are
    # int bitmasks restricted to P.
    pmask = 0
    for v in P:
        pmask |= 1 << v
    val = [0 * one, one]
    vs = [0, 0]
    dead = [True, False]
    for n in c.nodes[2:]:
        m = 0
        if n[0] == 'and':
            unit_set = set(n[2])
            is_dead = any(dead[k] for k in n[1]) or any(-l in unit_set for l in unit_set)
            x = 0 * one
            if not is_dead:
                x = one
                for k in n[1]:
                    x *= val[k]
                    m |= vs[k]
                for l in n[2]:
                    x *= w(l)
                    m |= (1 << abs(l)) & pmask
                for v in n[3]:
                    if v in P:
                        x *= both(v)
                        m |= 1 << v
        else:
            arcs = []
            for lits, ch in n[1]:
                ls = set(lits)
                if dead[ch] or any(-l in ls for l in ls):
                    continue
                am = vs[ch]
                for l in lits:
                    am |= (1 << abs(l)) & pmask
                m |= am
                arcs.append((lits, ch, am))
            is_dead = not arcs
            x = 0 * one
            for lits, ch, am in arcs:
                y = val[ch]
                for l in lits:
                    y *= w(l)
                miss = m & ~am
                while miss:
                    low = miss & -miss
                    y *= both(low.bit_length() - 1)
                    miss ^= low
                x += y
        if is_dead:
            m = 0
        val.append(x)
        vs.append(m)
        dead.append(is_dead)
    r = val[c.root]
    if not dead[c.root]:
        for v in P:
            if not vs[c.root] >> v & 1:
                r *= both(v)
    return r


def to_d4_model(c):
    """(nodes, arcs, root) in the d4-shaped dict model of ddnnf_verify.py: free
    var v -> shared synthetic node -v = OR(v -> T, -v -> T), units -> one arc to T."""
    nodes = {0: 'f', 1: 't'}
    arcs = {0: [], 1: []}
    for nid in range(2, len(c.nodes)):
        n = c.nodes[nid]
        if n[0] == 'and':
            nodes[nid] = 'a'
            a = [(k, []) for k in n[1]]
            if n[2]:
                a.append((1, list(n[2])))
            for v in n[3]:
                nodes[-v] = 'o'
                arcs[-v] = [(1, [v]), (1, [-v])]
                a.append((-v, []))
            arcs[nid] = a
        else:
            nodes[nid] = 'o'
            arcs[nid] = [(ch, list(l)) for l, ch in n[1]]
    return nodes, arcs, c.root


def dump(c, out):
    def ref(r):
        return 'F' if r == 0 else 'T' if r == 1 else str(r)
    out.write("bdnnf 1 nvars %d\n" % c.nvars)
    for k, v in c.header_kv:
        out.write("meta %s %s\n" % (k, v))
    for nid in range(2, len(c.nodes)):
        n = c.nodes[nid]
        if n[0] == 'and':
            s = "%d and" % nid + "".join(" " + ref(k) for k in n[1])
            if n[2]:
                s += " units" + "".join(" %d" % l for l in n[2])
            if n[3]:
                s += " free" + "".join(" %d" % v for v in n[3])
        else:
            s = "%d or" % nid + " ;".join(
                "".join(" %d" % l for l in lits) + " -> " + ref(ch) for lits, ch in n[1])
        out.write(s + "\n")
    out.write("root %s\n" % ref(c.root))
    out.write("props" + "".join(" " + p for i, p in enumerate(PROPS) if c.props >> i & 1) + "\n")
    for k, v in c.footer_kv:
        out.write("meta %s %s\n" % (k, v))


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["check", "count", "dump"])
    ap.add_argument("file")
    ap.add_argument("--cnf", help="CNF for weights (c p weight) and projection (c p show)")
    ap.add_argument("--structural", action="store_true",
                    help="print the structural count (literals ignored), which needs no "
                         "decomposability, instead of the weighted count")
    args = ap.parse_args()
    try:
        c = read(args.file)
    except FormatError as e:
        print("bdnnf: FORMAT ERROR: %s" % e)
        return 1
    if args.cmd == "check":
        errs = check_props(c)
        for e in errs:
            print("bdnnf: PROPERTY ERROR: %s" % e)
        if errs:
            return 1
        codecs = sorted({fr[4] for fr in c.frames})
        print("bdnnf: OK nodes=%d edges=%d lits=%d frames=%d codecs=%s props=%s"
              % (len(c.nodes) - 2, c.nedges, c.nlits, len(c.frames), codecs,
                 "|".join(p for i, p in enumerate(PROPS) if c.props >> i & 1) or "-"))
    elif args.cmd == "count" and args.structural:
        print(structural_count(c))
    elif args.cmd == "count":
        weights = proj = None
        if args.cnf:
            _, weights, proj = read_cnf_weights(args.cnf)
        print(weighted_count(c, weights, proj))
    else:
        dump(c, sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main())
