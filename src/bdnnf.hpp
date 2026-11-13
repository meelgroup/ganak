/******************************************
Copyright (C) 2023 Authors of GANAK, see AUTHORS file

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

#pragma once

// Header-only reader/writer for the BDNNF binary d-DNNF format (spec-bdnnf.md).
// Define USE_ZSTD and link libzstd for compressed frames.

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <iterator>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifdef USE_ZSTD
#include <zstd.h>
#endif

namespace bdnnf {

constexpr uint8_t MAGIC[8] = {0x89, 'B', 'D', 'N', 'N', 'F', '\r', '\n'};
constexpr uint8_t END_MAGIC[8] = {'B', 'D', 'N', 'N', 'F', 'E', 'N', 'D'};
constexpr uint32_t VERSION = 1;
constexpr size_t FRAME_HEAD = 20;
constexpr size_t FRAME_TAIL = 4;
constexpr size_t TRAILER = 16;

enum : uint8_t { FRAME_NODES = 1, FRAME_FOOTER = 2, FRAME_DEFS = 3 };
enum : uint8_t { CODEC_RAW = 0, CODEC_ZSTD = 1 };
enum : uint8_t { K_AND = 0, K_OR = 1 };
enum : uint8_t { TAG_KIND = 3, TAG_UNITS = 4, TAG_FREE = 8 };
enum : uint64_t {
  P_DECOMPOSABLE = 1, P_DETERMINISTIC = 2, P_SMOOTH = 4, P_CLEAN = 8, P_ALL = 15
};

struct Error : std::runtime_error {
  using std::runtime_error::runtime_error;
};

using KV = std::vector<std::pair<std::string, std::string>>;

inline const uint32_t* crc_tables() {
  static const std::vector<uint32_t> t = [] {
    std::vector<uint32_t> tab(8 * 256);
    for (uint32_t i = 0; i < 256; i++) {
      uint32_t c = i;
      for (int k = 0; k < 8; k++) c = (c >> 1) ^ (0x82F63B78u & (0u - (c & 1u)));
      tab[i] = c;
    }
    for (int s = 1; s < 8; s++)
      for (uint32_t i = 0; i < 256; i++) {
        uint32_t prev = tab[(s - 1) * 256 + i];
        tab[s * 256 + i] = (prev >> 8) ^ tab[prev & 0xff];
      }
    return tab;
  }();
  return t.data();
}

// CRC-32C, slicing-by-8.
inline uint32_t crc32c(const uint8_t* p, size_t n) {
  const uint32_t* t = crc_tables();
  uint32_t crc = ~0u;
  auto ld = [](const uint8_t* q) {
    return (uint32_t)q[0] | ((uint32_t)q[1] << 8) | ((uint32_t)q[2] << 16) | ((uint32_t)q[3] << 24);
  };
  while (n >= 8) {
    uint32_t lo = ld(p) ^ crc;
    uint32_t hi = ld(p + 4);
    crc = t[7 * 256 + (lo & 0xff)] ^ t[6 * 256 + ((lo >> 8) & 0xff)] ^
          t[5 * 256 + ((lo >> 16) & 0xff)] ^ t[4 * 256 + (lo >> 24)] ^
          t[3 * 256 + (hi & 0xff)] ^ t[2 * 256 + ((hi >> 8) & 0xff)] ^
          t[1 * 256 + ((hi >> 16) & 0xff)] ^ t[0 * 256 + (hi >> 24)];
    p += 8;
    n -= 8;
  }
  while (n--) crc = (crc >> 8) ^ t[(crc ^ *p++) & 0xff];
  return ~crc;
}

inline void put_varint(std::vector<uint8_t>& b, uint64_t v) {
  while (v >= 0x80) {
    b.push_back((uint8_t)(v | 0x80));
    v >>= 7;
  }
  b.push_back((uint8_t)v);
}
inline void put_le(std::vector<uint8_t>& b, uint64_t v, int bytes) {
  for (int i = 0; i < bytes; i++) b.push_back((uint8_t)(v >> (8 * i)));
}
inline uint64_t get_le(const uint8_t* p, int bytes) {
  uint64_t v = 0;
  for (int i = 0; i < bytes; i++) v |= (uint64_t)p[i] << (8 * i);
  return v;
}
inline void put_str(std::vector<uint8_t>& b, const std::string& s) {
  put_varint(b, s.size());
  b.insert(b.end(), s.begin(), s.end());
}
inline void put_kv(std::vector<uint8_t>& b, const KV& kv) {
  put_varint(b, kv.size());
  for (const auto& [k, v] : kv) {
    put_str(b, k);
    put_str(b, v);
  }
}

struct Cursor {
  const uint8_t* p;
  const uint8_t* end;
  Cursor(const uint8_t* b, size_t n) : p(b), end(b + n) {}
  [[nodiscard]] bool done() const { return p == end; }
  uint8_t byte() {
    if (p == end) throw Error("truncated record");
    return *p++;
  }
  uint64_t varint() {
    uint64_t v = 0;
    for (int sh = 0;; sh += 7) {
      if (p == end) throw Error("truncated varint");
      const uint8_t c = *p++;
      if (sh == 63 && c > 1) throw Error("varint overflows 64 bits");
      v |= (uint64_t)(c & 0x7f) << sh;
      if (!(c & 0x80)) return v;
      if (sh >= 63) throw Error("varint too long");
    }
  }
  std::string str() {
    const uint64_t n = varint();
    if (n > (uint64_t)(end - p)) throw Error("truncated string");
    std::string s((const char*)p, (size_t)n);
    p += n;
    return s;
  }
  KV kv() {
    const uint64_t n = varint();
    KV r;
    for (uint64_t i = 0; i < n; i++) {
      std::string k = str();
      std::string v = str();
      r.emplace_back(std::move(k), std::move(v));
    }
    return r;
  }
};

inline uint64_t lit_code(int32_t l) {
  const uint64_t v = l > 0 ? (uint64_t)l : (uint64_t)(-(int64_t)l);
  return 2 * v + (l < 0 ? 1 : 0);
}
inline int32_t code_lit(uint64_t c) {
  const int64_t v = (int64_t)(c >> 1);
  return (int32_t)((c & 1) ? -v : v);
}

struct Arc {
  std::vector<int32_t> lits;
  uint64_t child = 0;
};

// The definitions frame: AIG functions of vars not in the circuit. A ref r is
// k = r >> 1, negated if r & 1: k = 0 is TRUE, 1..nvars a var, nvars + 1 + j
// the j-th AND node.
struct Defs {
  std::vector<std::pair<uint64_t, uint64_t>> nodes;  // AND(ref, ref)
  std::vector<std::pair<uint32_t, uint64_t>> defs;   // (var, ref)
  [[nodiscard]] bool empty() const { return defs.empty(); }
};

inline std::vector<uint8_t> encode_defs(const Defs& d) {
  std::vector<uint8_t> p;
  put_varint(p, d.nodes.size());
  for (const auto& [a, b] : d.nodes) { put_varint(p, a); put_varint(p, b); }
  put_varint(p, d.defs.size());
  for (const auto& [v, r] : d.defs) { put_varint(p, v); put_varint(p, r); }
  return p;
}

// Also checks: refs in range, AND nodes in topological order, each var defined
// once. (Acyclicity through the defined vars is left to the evaluator.)
inline void decode_defs(const std::vector<uint8_t>& raw, uint32_t nvars, Defs& d) {
  Cursor c(raw.data(), raw.size());
  d.nodes.clear();
  d.defs.clear();
  const uint64_t nn = c.varint();
  if (nn > raw.size()) throw Error("bad definitions frame");
  for (uint64_t j = 0; j < nn; j++) {
    const uint64_t a = c.varint(), b = c.varint();
    if ((a >> 1) >= (uint64_t)nvars + 1 + j || (b >> 1) >= (uint64_t)nvars + 1 + j)
      throw Error("definition AIG node refers to a later node");
    d.nodes.emplace_back(a, b);
  }
  const uint64_t nd = c.varint();
  std::vector<char> seen((size_t)nvars + 1, 0);
  for (uint64_t i = 0; i < nd; i++) {
    const uint64_t v = c.varint(), r = c.varint();
    if (v < 1 || v > nvars) throw Error("defined var out of range");
    if (seen[v]) throw Error("var defined twice");
    seen[v] = 1;
    if ((r >> 1) > (uint64_t)nvars + nn) throw Error("definition ref out of range");
    d.defs.emplace_back((uint32_t)v, r);
  }
  if (!c.done()) throw Error("trailing bytes in definitions frame");
}

// Sets every defined var in val (by var: 1 true, -1 false, 0 unset) from the
// others. Undefined vars left unset count as false. Throws on a cycle.
inline void eval_defs(const Defs& d, uint32_t nvars, std::vector<signed char>& val) {
  std::vector<uint64_t> def_of((size_t)nvars + 1, UINT64_MAX);
  for (const auto& [v, r] : d.defs) def_of[v] = r;
  std::vector<signed char> node_val(d.nodes.size(), 0);
  std::vector<char> state((size_t)nvars + 1, 0);  // 1: being evaluated, 2: done
  // Iterative over refs: a work stack of (ref, expanded).
  auto value = [&](uint64_t ref0) -> bool {
    std::vector<std::pair<uint64_t, bool>> st{{ref0 & ~1ULL, false}};
    while (!st.empty()) {
      auto [r, expanded] = st.back();
      const uint64_t k = r >> 1;
      if (k == 0) { st.pop_back(); continue; }
      if (k <= nvars) {
        const uint32_t v = (uint32_t)k;
        if (def_of[v] == UINT64_MAX || state[v] == 2) { st.pop_back(); continue; }
        if (!expanded) {
          if (state[v] == 1) throw Error("definitions are cyclic");
          state[v] = 1;
          st.back().second = true;
          st.push_back({def_of[v] & ~1ULL, false});
          continue;
        }
        st.pop_back();
        const uint64_t dr = def_of[v];
        const uint64_t dk = dr >> 1;
        bool x = dk == 0 ? true : dk <= nvars ? val[dk] > 0 : node_val[dk - nvars - 1] > 0;
        if (dr & 1) x = !x;
        val[v] = x ? 1 : -1;
        state[v] = 2;
        continue;
      }
      const size_t j = (size_t)(k - nvars - 1);
      if (node_val[j]) { st.pop_back(); continue; }
      const auto [a, b] = d.nodes[j];
      if (!expanded) {
        st.back().second = true;
        st.push_back({a & ~1ULL, false});
        st.push_back({b & ~1ULL, false});
        continue;
      }
      st.pop_back();
      auto lit = [&](uint64_t x) {
        const uint64_t xk = x >> 1;
        bool y = xk == 0 ? true : xk <= nvars ? val[xk] > 0 : node_val[xk - nvars - 1] > 0;
        return (x & 1) ? !y : y;
      };
      node_val[j] = (lit(a) && lit(b)) ? 1 : -1;
    }
    const uint64_t k0 = ref0 >> 1;
    bool y = k0 == 0 ? true : k0 <= nvars ? val[k0] > 0 : node_val[k0 - nvars - 1] > 0;
    return (ref0 & 1) ? !y : y;
  };
  for (const auto& [v, r] : d.defs) {
    if (state[v] == 2) continue;
    (void)value(((uint64_t)v) << 1);
  }
}

// Decoded node. Vectors are reused across decode() calls to avoid allocation.
struct Node {
  uint8_t kind = K_AND;
  std::vector<uint64_t> kids;
  std::vector<int32_t> units;
  std::vector<uint32_t> free;
  std::vector<Arc> arcs;
};

class Writer {
public:
  int zstd_level = 1;  // 0 = always raw
  size_t frame_target = 128 * 1024;

  Writer() = default;
  Writer(const Writer&) = delete;
  Writer& operator=(const Writer&) = delete;
  ~Writer() {
#ifdef USE_ZSTD
    if (cctx) ZSTD_freeCCtx(cctx);
#endif
    if (f && f != stdout) fclose(f);
  }

  void open(const std::string& path, uint32_t _nvars, const KV& kv = {}) {
    nvars = _nvars;
    if (path == "-") f = stdout;
    else f = fopen(path.c_str(), "wb");
    if (!f) throw Error("cannot open for writing: " + path);
    setvbuf(f, nullptr, _IOFBF, 1 << 20);
    std::vector<uint8_t> payload;
    put_varint(payload, nvars);
    put_kv(payload, kv);
    std::vector<uint8_t> h(MAGIC, MAGIC + 8);
    put_le(h, VERSION, 4);
    put_le(h, payload.size(), 4);
    h.insert(h.end(), payload.begin(), payload.end());
    write_bytes(h.data(), h.size());
  }

  [[nodiscard]] uint64_t next_id() const { return next; }
  [[nodiscard]] uint64_t num_nodes() const { return next - 2; }
  [[nodiscard]] uint64_t num_edges() const { return nedges; }

  template <class Kids, class Units, class Free>
  uint64_t add_and(const Kids& kids, const Units& units = std::vector<int>{},
                   const Free& free = std::vector<int>{}) {
    const uint64_t self = next;
    uint8_t tag = K_AND;
    if (!std::empty(units)) tag |= TAG_UNITS;
    if (!std::empty(free)) tag |= TAG_FREE;
    buf.push_back(tag);
    put_varint(buf, std::size(kids));
    for (auto c : kids) put_ref(self, (uint64_t)c);
    if (tag & TAG_UNITS) put_lits(units);
    if (tag & TAG_FREE) put_vars(free);
    nedges += std::size(kids);
    return finish_node();
  }

  void begin_or() {
    or_arcs = 0;
    or_buf.clear();
  }
  template <class Lits>
  void or_arc(const Lits& lits, uint64_t child) {
    const uint64_t self = next;
    std::swap(buf, or_buf);
    put_lits(lits);
    put_ref(self, child);
    std::swap(buf, or_buf);
    or_arcs++;
  }
  uint64_t end_or() {
    buf.push_back(K_OR);
    put_varint(buf, or_arcs);
    buf.insert(buf.end(), or_buf.begin(), or_buf.end());
    nedges += or_arcs;
    return finish_node();
  }

  void finish(uint64_t root, uint64_t props, const KV& kv = {}, const Defs* defs = nullptr) {
    if (root >= next) throw Error("root id out of range");
    if (props & ~(uint64_t)P_ALL) throw Error("unknown property bits");
    if (frame_nodes) flush_frame();
    uint64_t defs_size = 0;
    if (defs && !defs->empty()) {
      const auto d = encode_defs(*defs);
      const uint64_t before = offset;
      write_frame(FRAME_DEFS, CODEC_RAW, 0, d.data(), d.size(), d.size(), crc32c(d.data(), d.size()));
      defs_size = offset - before;
    }
    std::vector<uint8_t> p;
    put_varint(p, next - 2);
    put_varint(p, nedges);
    put_varint(p, nlits);
    put_varint(p, root);
    put_varint(p, props);
    put_varint(p, index.size());
    for (const auto& [sz, nn] : index) {
      put_varint(p, sz);
      put_varint(p, nn);
    }
    put_varint(p, defs_size);
    put_kv(p, kv);
    const uint64_t footer_off = offset;
    write_frame(FRAME_FOOTER, CODEC_RAW, 0, p.data(), p.size(), p.size(), crc32c(p.data(), p.size()));
    std::vector<uint8_t> t;
    put_le(t, footer_off, 8);
    t.insert(t.end(), END_MAGIC, END_MAGIC + 8);
    write_bytes(t.data(), t.size());
    if (fflush(f) != 0) throw Error("write error");
    if (f != stdout) {
      const int r = fclose(f);
      f = nullptr;
      if (r != 0) throw Error("write error on close");
    }
  }

private:
  FILE* f = nullptr;
  uint32_t nvars = 0;
  uint64_t next = 2;
  uint64_t offset = 0;
  uint64_t nedges = 0;
  uint64_t nlits = 0;
  uint32_t frame_nodes = 0;
  std::vector<uint8_t> buf;
  std::vector<uint8_t> or_buf;
  uint64_t or_arcs = 0;
  std::vector<uint64_t> codes;
  std::vector<std::pair<uint64_t, uint64_t>> index;  // (frame size, nnodes)
  std::vector<uint8_t> zbuf;
#ifdef USE_ZSTD
  ZSTD_CCtx* cctx = nullptr;
#endif

  void write_bytes(const uint8_t* p, size_t n) {
    if (n && fwrite(p, 1, n, f) != n) throw Error("write error");
    offset += n;
  }

  void put_ref(uint64_t self, uint64_t child) {
    if (child <= 1) put_varint(buf, child);
    else if (child >= self) throw Error("child id not smaller than parent");
    else put_varint(buf, self - child + 1);
  }

  void put_codes() {
    std::sort(codes.begin(), codes.end());
    codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
    put_varint(buf, codes.size());
    for (size_t i = 0; i < codes.size(); i++)
      put_varint(buf, i == 0 ? codes[0] : codes[i] - codes[i - 1] - 1);
    nlits += codes.size();
  }
  template <class Lits>
  void put_lits(const Lits& lits) {
    codes.clear();
    for (auto l : lits) {
      const int64_t v = l > 0 ? (int64_t)l : -(int64_t)l;
      if (l == 0 || v > (int64_t)nvars) throw Error("literal out of range");
      codes.push_back(lit_code((int32_t)l));
    }
    put_codes();
  }
  template <class Vars>
  void put_vars(const Vars& vars) {
    codes.clear();
    for (auto v : vars) {
      if ((int64_t)v <= 0 || (uint64_t)v > nvars) throw Error("free variable out of range");
      codes.push_back((uint64_t)v);
    }
    put_codes();
  }

  uint64_t finish_node() {
    frame_nodes++;
    if (buf.size() >= frame_target || frame_nodes == UINT32_MAX) flush_frame();
    return next++;
  }

  void flush_frame() {
    const uint32_t crc = crc32c(buf.data(), buf.size());
    if (buf.size() > UINT32_MAX) throw Error("frame too large");
#ifdef USE_ZSTD
    if (zstd_level > 0) {
      if (!cctx) cctx = ZSTD_createCCtx();
      zbuf.resize(ZSTD_compressBound(buf.size()));
      const size_t r = ZSTD_compressCCtx(cctx, zbuf.data(), zbuf.size(), buf.data(), buf.size(), zstd_level);
      if (ZSTD_isError(r)) throw Error(std::string("zstd: ") + ZSTD_getErrorName(r));
      if (r < buf.size()) {
        write_frame(FRAME_NODES, CODEC_ZSTD, frame_nodes, zbuf.data(), r, buf.size(), crc);
        reset_frame();
        return;
      }
    }
#endif
    write_frame(FRAME_NODES, CODEC_RAW, frame_nodes, buf.data(), buf.size(), buf.size(), crc);
    reset_frame();
  }
  void reset_frame() {
    buf.clear();
    frame_nodes = 0;
  }

  void write_frame(uint8_t kind, uint8_t codec, uint32_t nn, const uint8_t* data,
                   size_t stored, size_t raw, uint32_t crc) {
    std::vector<uint8_t> h;
    h.push_back(kind);
    h.push_back(codec);
    put_le(h, 0, 2);
    put_le(h, nn, 4);
    put_le(h, stored, 4);
    put_le(h, raw, 4);
    put_le(h, crc, 4);
    write_bytes(h.data(), h.size());
    write_bytes(data, stored);
    h.clear();
    put_le(h, stored, 4);
    write_bytes(h.data(), h.size());
    if (kind == FRAME_NODES) index.emplace_back(stored + FRAME_HEAD + FRAME_TAIL, nn);
  }
};

struct FrameInfo {
  uint64_t offset;
  uint64_t size;
  uint64_t first_id;
  uint32_t nnodes;
};

// Decodes one node record at `c`, for node id `self`.
inline void decode_node(Cursor& c, uint64_t self, uint32_t nvars, Node& n) {
  const uint8_t tag = c.byte();
  n.kind = tag & TAG_KIND;
  n.kids.clear();
  n.units.clear();
  n.free.clear();
  n.arcs.clear();
  auto ref = [&]() -> uint64_t {
    const uint64_t r = c.varint();
    if (r <= 1) return r;
    if (r > self - 1) throw Error("child reference out of range at node " + std::to_string(self));
    return self - r + 1;
  };
  auto codes = [&](auto&& sink) {
    const uint64_t k = c.varint();
    uint64_t prev = 0;
    for (uint64_t i = 0; i < k; i++) {
      const uint64_t g = c.varint();
      const uint64_t cur = i == 0 ? g : prev + g + 1;
      if (i && cur <= prev) throw Error("code overflow");
      sink(cur);
      prev = cur;
    }
  };
  auto lits = [&](std::vector<int32_t>& out) {
    codes([&](uint64_t code) {
      if (code < 2 || (code >> 1) > nvars) throw Error("literal out of range at node " + std::to_string(self));
      out.push_back(code_lit(code));
    });
  };
  if (n.kind == K_AND) {
    if (tag & ~(TAG_KIND | TAG_UNITS | TAG_FREE)) throw Error("reserved tag bits set");
    const uint64_t nk = c.varint();
    if (nk > (uint64_t)(c.end - c.p)) throw Error("truncated AND");
    for (uint64_t i = 0; i < nk; i++) n.kids.push_back(ref());
    if (tag & TAG_UNITS) lits(n.units);
    if (tag & TAG_FREE)
      codes([&](uint64_t v) {
        if (v < 1 || v > nvars) throw Error("free variable out of range at node " + std::to_string(self));
        n.free.push_back((uint32_t)v);
      });
  } else if (n.kind == K_OR) {
    if (tag & ~TAG_KIND) throw Error("reserved tag bits set");
    const uint64_t na = c.varint();
    if (na > (uint64_t)(c.end - c.p)) throw Error("truncated OR");
    n.arcs.resize(na);
    for (uint64_t i = 0; i < na; i++) {
      n.arcs[i].lits.clear();
      lits(n.arcs[i].lits);
      n.arcs[i].child = ref();
    }
  } else {
    throw Error("unknown node kind at node " + std::to_string(self));
  }
}

class Reader {
public:
  uint32_t nvars = 0;
  KV header_kv;
  uint64_t nnodes = 0, nedges = 0, nlits = 0, root = 0, props = 0;
  KV footer_kv;
  std::vector<FrameInfo> frames;
  bool have_footer = false;
  Defs defs;

  Reader() = default;
  Reader(const Reader&) = delete;
  Reader& operator=(const Reader&) = delete;
  ~Reader() {
    if (f && f != stdin) fclose(f);
  }

  // "-" reads stdin as a stream: only scan() works, the footer arrives at its end.
  void open(const std::string& path) {
    if (path == "-") {
      f = stdin;
      seekable = false;
    } else {
      f = fopen(path.c_str(), "rb");
      if (!f) throw Error("cannot open: " + path);
      seekable = true;
    }
    read_header();
    if (seekable) read_footer_via_trailer();
  }

  [[nodiscard]] std::string get(const KV& kv, const std::string& key) const {
    for (const auto& [k, v] : kv) if (k == key) return v;
    return "";
  }

  // Calls fn(id, node) for every node in id order, then reads the footer.
  template <class Fn>
  void scan(Fn&& fn) {
    defs_seen = false;
    if (seekable) seek(data_start);
    pos = data_start;
    uint64_t id = 2;
    size_t fi = 0;
    while (true) {
      uint8_t kind;
      uint32_t nn;
      const uint64_t off = pos;
      read_frame_at_pos(kind, nn, raw);
      if (kind == FRAME_DEFS) {
        if (have_footer && off + defs_size != footer_off) throw Error("definitions frame not right before the footer");
        decode_defs(raw, nvars, defs);
        defs_seen = true;
        continue;
      }
      if (kind == FRAME_FOOTER) {
        if (have_footer) {
          if (off != footer_off) throw Error("footer offset mismatch");
        } else {
          parse_footer(raw, off);
          read_trailer_stream(off);
        }
        break;
      }
      if (have_footer && (fi >= frames.size() || frames[fi].offset != off || frames[fi].nnodes != nn))
        throw Error("frame does not match footer index");
      fi++;
      Cursor c(raw.data(), raw.size());
      for (uint32_t i = 0; i < nn; i++, id++) {
        decode_node(c, id, nvars, node);
        fn(id, (const Node&)node);
      }
      if (!c.done()) throw Error("trailing bytes in frame");
    }
    if (id != nnodes + 2 || fi != frames.size()) throw Error("node count does not match footer");
    if ((defs_size != 0) != defs_seen) throw Error("definitions frame does not match the footer");
  }

  // Calls fn(id, node) for every node in reverse id order. Seekable files only.
  template <class Fn>
  void scan_backward(Fn&& fn) {
    if (!seekable) throw Error("backward scan needs a seekable file");
    std::vector<size_t> starts;
    for (size_t fi = frames.size(); fi-- > 0;) {
      load_frame(fi, raw);
      starts.clear();
      Cursor c(raw.data(), raw.size());
      for (uint32_t i = 0; i < frames[fi].nnodes; i++) {
        starts.push_back((size_t)(c.p - raw.data()));
        decode_node(c, frames[fi].first_id + i, nvars, node);
      }
      for (size_t i = starts.size(); i-- > 0;) {
        Cursor ci(raw.data() + starts[i], raw.size() - starts[i]);
        decode_node(ci, frames[fi].first_id + i, nvars, node);
        fn(frames[fi].first_id + i, (const Node&)node);
      }
    }
  }

  // Index of the frame holding node `id`.
  [[nodiscard]] size_t frame_of(uint64_t id) const {
    if (id < 2 || id >= nnodes + 2) throw Error("node id out of range");
    auto it = std::upper_bound(frames.begin(), frames.end(), id,
        [](uint64_t x, const FrameInfo& fr) { return x < fr.first_id; });
    return (size_t)(it - frames.begin()) - 1;
  }

  void load_frame(size_t fi, std::vector<uint8_t>& out) {
    if (!seekable) throw Error("random access needs a seekable file");
    seek(frames[fi].offset);
    pos = frames[fi].offset;
    uint8_t kind;
    uint32_t nn;
    read_frame_at_pos(kind, nn, out);
    if (kind != FRAME_NODES || nn != frames[fi].nnodes) throw Error("frame does not match footer index");
  }

private:
  FILE* f = nullptr;
  bool seekable = false;
  uint64_t pos = 0;
  uint64_t data_start = 0;
  uint64_t footer_off = 0;
  uint64_t defs_size = 0;
  bool defs_seen = false;
  std::vector<uint8_t> raw, stored;
  Node node;

  void read_exact(uint8_t* p, size_t n) {
    if (n && fread(p, 1, n, f) != n) throw Error("unexpected end of file");
    pos += n;
  }
  void seek(uint64_t off) {
    if (fseeko(f, (off_t)off, SEEK_SET) != 0) throw Error("seek failed");
  }

  void read_header() {
    uint8_t h[16];
    read_exact(h, 16);
    if (memcmp(h, MAGIC, 8) != 0) throw Error("not a BDNNF file (bad magic)");
    const uint32_t ver = (uint32_t)get_le(h + 8, 4);
    if (ver != VERSION) throw Error("unsupported BDNNF version " + std::to_string(ver));
    const uint32_t hlen = (uint32_t)get_le(h + 12, 4);
    std::vector<uint8_t> p(hlen);
    read_exact(p.data(), hlen);
    Cursor c(p.data(), p.size());
    const uint64_t nv = c.varint();
    if (nv > (uint64_t)INT32_MAX) throw Error("nvars too large");
    nvars = (uint32_t)nv;
    header_kv = c.kv();
    if (!c.done()) throw Error("trailing bytes in header");
    data_start = pos;
  }

  void read_frame_at_pos(uint8_t& kind, uint32_t& nn, std::vector<uint8_t>& out) {
    uint8_t h[FRAME_HEAD];
    read_exact(h, FRAME_HEAD);
    kind = h[0];
    const uint8_t codec = h[1];
    if (kind != FRAME_NODES && kind != FRAME_FOOTER && kind != FRAME_DEFS) throw Error("bad frame kind");
    if (get_le(h + 2, 2) != 0) throw Error("reserved frame bits set");
    nn = (uint32_t)get_le(h + 4, 4);
    const uint32_t slen = (uint32_t)get_le(h + 8, 4);
    const uint32_t rlen = (uint32_t)get_le(h + 12, 4);
    const uint32_t crc = (uint32_t)get_le(h + 16, 4);
    if (codec == CODEC_RAW) {
      if (slen != rlen) throw Error("raw frame length mismatch");
      out.resize(rlen);
      read_exact(out.data(), rlen);
    } else if (codec == CODEC_ZSTD) {
#ifdef USE_ZSTD
      stored.resize(slen);
      read_exact(stored.data(), slen);
      out.resize(rlen);
      const size_t r = ZSTD_decompress(out.data(), rlen, stored.data(), slen);
      if (ZSTD_isError(r) || r != rlen) throw Error("zstd frame failed to decompress");
#else
      throw Error("file has zstd frames, but this tool was built without zstd");
#endif
    } else {
      throw Error("unknown codec " + std::to_string(codec));
    }
    uint8_t t[FRAME_TAIL];
    read_exact(t, FRAME_TAIL);
    if (get_le(t, 4) != slen) throw Error("frame tail mismatch");
    if (crc32c(out.data(), out.size()) != crc) throw Error("frame CRC mismatch");
  }

  void parse_footer(const std::vector<uint8_t>& p, uint64_t off) {
    Cursor c(p.data(), p.size());
    nnodes = c.varint();
    nedges = c.varint();
    nlits = c.varint();
    root = c.varint();
    props = c.varint();
    if (props & ~(uint64_t)P_ALL) throw Error("unknown property bits");
    if (root >= nnodes + 2) throw Error("root out of range");
    const uint64_t nf = c.varint();
    frames.clear();
    uint64_t o = data_start, id = 2;
    for (uint64_t i = 0; i < nf; i++) {
      const uint64_t sz = c.varint();
      const uint64_t nn = c.varint();
      if (nn == 0 || nn > UINT32_MAX || sz < FRAME_HEAD + FRAME_TAIL) throw Error("bad frame index");
      frames.push_back({o, sz, id, (uint32_t)nn});
      o += sz;
      id += nn;
    }
    defs_size = c.varint();
    if (o + defs_size != off) throw Error("frame index does not add up to footer offset");
    if (id != nnodes + 2) throw Error("frame index does not add up to node count");
    footer_kv = c.kv();
    if (!c.done()) throw Error("trailing bytes in footer");
    footer_off = off;
    have_footer = true;
  }

  void read_trailer_stream(uint64_t off) {
    uint8_t t[TRAILER];
    read_exact(t, TRAILER);
    if (memcmp(t + 8, END_MAGIC, 8) != 0 || get_le(t, 8) != off) throw Error("bad trailer");
  }

  void read_footer_via_trailer() {
    if (fseeko(f, 0, SEEK_END) != 0) throw Error("seek failed");
    const off_t sz = ftello(f);
    if (sz < (off_t)(data_start + FRAME_HEAD + FRAME_TAIL + TRAILER))
      throw Error("file too short: no footer (incomplete file?)");
    seek((uint64_t)sz - TRAILER);
    uint8_t t[TRAILER];
    read_exact(t, TRAILER);
    if (memcmp(t + 8, END_MAGIC, 8) != 0) throw Error("no trailer: file is incomplete");
    const uint64_t off = get_le(t, 8);
    if (off < data_start || off + FRAME_HEAD + FRAME_TAIL + TRAILER > (uint64_t)sz)
      throw Error("bad footer offset");
    seek(off);
    pos = off;
    uint8_t kind;
    uint32_t nn;
    read_frame_at_pos(kind, nn, raw);
    if (kind != FRAME_FOOTER) throw Error("trailer does not point at the footer");
    if (pos + TRAILER != (uint64_t)sz) throw Error("garbage between footer and trailer");
    parse_footer(raw, off);
    if (defs_size) {
      seek(off - defs_size);
      pos = off - defs_size;
      read_frame_at_pos(kind, nn, raw);
      if (kind != FRAME_DEFS || pos != off) throw Error("bad definitions frame");
      decode_defs(raw, nvars, defs);
    }
  }
};

}  // namespace bdnnf
