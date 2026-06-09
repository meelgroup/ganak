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

// Observer the counting engine notifies at every search event relevant to
// d-DNNF compilation. Default NullCompiler is a no-op; the public hooks are
// non-virtual inline wrappers that short-circuit on `active_`, so with
// --compile off the steady-state loop pays only one not-taken branch per
// event. With --compile on, DDNNFCompiler (compiler.cpp) drives a DDNNFCircuit.
struct Compiler {
  virtual ~Compiler() = default;

  [[nodiscard]] bool active() const { return active_; }
  // Compiler needs the full DPLL tree, so it disallows the small-formula CMS shortcut.
  [[nodiscard]] bool allows_cms_shortcut() const { return !active_; }

  // ---- search events ----
  void new_decision_level()    { if (active_) do_new_decision_level(); }
  void save_left_lits()        { if (active_) do_save_left_lits(); }
  void build_level_node()      { if (active_) do_build_level_node(); }
  void record_comp_node()      { if (active_) do_record_comp_node(); }
  void attach_to_parent()      { if (active_) do_attach_to_parent(); }
  void branch_reset()          { if (active_) do_branch_reset(); }
  void free_var(uint32_t v)    { if (active_) do_free_var(v); }
  void cache_hit(int node)     { if (active_) do_cache_hit(node); }
  void sat_witness_capture(int sat_start_dec_level) {
    if (active_) do_sat_witness_capture(sat_start_dec_level);
  }
  void sat_witness_apply()     { if (active_) do_sat_witness_apply(); }

  // ---- lifecycle ----
  void finalize_root()                 { if (active_) do_finalize_root(); }
  void finalize_and_write(bool unsat)  { if (active_) do_finalize_and_write(unsat); }

protected:
  // Set to true by DDNNFCompiler; stays false for NullCompiler.
  bool active_ = false;

  // Subclasses override these; only called when `active_` is true.
  virtual void do_new_decision_level() {}
  virtual void do_save_left_lits() {}
  virtual void do_build_level_node() {}
  virtual void do_record_comp_node() {}
  virtual void do_attach_to_parent() {}
  virtual void do_branch_reset() {}
  virtual void do_free_var(uint32_t /*v*/) {}
  virtual void do_cache_hit(int /*node*/) {}
  virtual void do_sat_witness_capture(int /*sat_start_dec_level*/) {}
  virtual void do_sat_witness_apply() {}
  virtual void do_finalize_root() {}
  virtual void do_finalize_and_write(bool /*unsat*/) {}
};

std::unique_ptr<Compiler> make_null_compiler();
std::unique_ptr<Compiler> make_ddnnf_compiler(Counter& counter, const std::string& fname);

}
