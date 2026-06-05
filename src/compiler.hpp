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

#include <cstdint>
#include <memory>
#include <string>

namespace GanakInt {

class Counter;

// Observer that the counting engine notifies at every search event relevant to
// d-DNNF compilation. The default (NullCompiler) does nothing, so the engine can
// call these hooks unconditionally -- no scattered `if (compiling())` checks. When
// `--compile` is on, a DDNNFCompiler (defined in compiler.cpp) drives a DDNNFCircuit
// from these events.
struct Compiler {
  virtual ~Compiler() = default;

  // ---- queries ----
  // Whether compilation is active (true only for the real compiler).
  virtual bool active() const { return false; }
  // Whether the engine may take the small-formula CMS shortcut (count_using_cms).
  // The compiler needs the full DPLL tree, so it disallows the shortcut.
  virtual bool allows_cms_shortcut() const { return true; }

  // ---- search events ----
  virtual void new_decision_level() {}      // a decision was made: open a fresh level
  virtual void save_left_lits() {}          // left branch done: snapshot its trail lits
  virtual void build_level_node() {}        // both branches done: build this level's OR node
  virtual void record_comp_node() {}        // store the built node against the cached component
  virtual void attach_to_parent() {}        // AND the built node into the parent branch
  virtual void branch_reset() {}            // current branch's count is discarded; drop its children
  virtual void free_var(uint32_t /*v*/) {}  // a free (unconstrained) var: factor of two
  virtual void cache_hit(int /*node*/) {}   // a cached component was reused: share its sub-DAG
  virtual void sat_witness_capture(int /*sat_start_dec_level*/) {} // SAT leaf: record synthesized-var witness
  virtual void sat_witness_apply() {}       // install the captured witness as this level's leaf

  // ---- lifecycle ----
  virtual void finalize_root() {}                  // EXIT state: build the root from level-0 children
  virtual void finalize_and_write(bool /*unsat*/) {} // write the d4 .nnf file
};

std::unique_ptr<Compiler> make_null_compiler();
std::unique_ptr<Compiler> make_ddnnf_compiler(Counter& counter, const std::string& fname);

}
