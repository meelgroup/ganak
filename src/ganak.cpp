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

#include "ganak.hpp"
#include "outer_counter.hpp"
#include <cstdlib>
#include <set>

using namespace GanakInt;

#if defined _WIN32
    #define DLL_PUBLIC __declspec(dllexport)
#else
    #define DLL_PUBLIC __attribute__ ((visibility ("default")))
    #define DLL_LOCAL  __attribute__ ((visibility ("hidden")))
#endif

struct CDat {
  CounterConfiguration conf;
  std::unique_ptr<CMSat::FieldGen> fg;
  set<uint32_t> indeps;
  set<uint32_t> opt_indeps;
  map<GanakInt::Lit, FF> lit_weights;
  uint32_t nvars = 0;
  vector<vector<GanakInt::Lit>> irred_cls;
  vector<uint8_t> irred_cls_used;
  vector<std::pair<vector<GanakInt::Lit>, uint32_t>> red_cls;
  bool print_indep_distrib = false;
  bool is_approximate = false;
  uint64_t max_num_cache_lookups = 0;
  uint64_t max_cache_elems = 0;
};

vector<vector<uint32_t>> find_disconnected(const CDat& dat);

DLL_PUBLIC Ganak::Ganak(CounterConfiguration& conf, FG& fg) :
    cdat(std::make_unique<CDat>()) {
  cdat->conf = conf;
  cdat->fg = fg->dup();
}

DLL_PUBLIC Ganak::~Ganak() = default;

DLL_PUBLIC FF Ganak::count(uint8_t bits_jobs, int num_threads) {
  auto cnt = cdat->fg->one();

  // Check for empty clause
  if (std::any_of(cdat->irred_cls.begin(), cdat->irred_cls.end(),
      [](const auto& cl) { return cl.empty(); })) {
    cout << "c o intermediate count: " << *cdat->fg->zero() << endl;
    return cdat->fg->zero();
  }

  auto bags = find_disconnected(*cdat);
  vector<int32_t> var_to_bag(cdat->nvars+1, -1);
  for(uint32_t i = 0; i < bags.size(); i++) {
    for(auto& v: bags[i]) {
      assert(v < cdat->nvars+1);
      assert(v > 0);
      assert(var_to_bag[v] == -1);
      var_to_bag[v] = i;
    }
  }
  vector<vector<vector<GanakInt::Lit>>> bag_to_irred_cls(bags.size());
  vector<vector<pair<vector<GanakInt::Lit>, uint32_t>>> bag_to_red_cls(bags.size());
  for(const auto& cl: cdat->irred_cls) {
    assert(!cl.empty() && "We filtered out empty clauses above");
    const int b = var_to_bag[cl[0].var()];
    bag_to_irred_cls[b].push_back(cl);
  }
  for(const auto& cl_p: cdat->red_cls) {
    assert(!cl_p.first.empty());
    const int b = var_to_bag[cl_p.first[0].var()];
    bool ok = true;
    for(const auto& l: cl_p.first) {
      if (var_to_bag[l.var()] != b) { ok = false; break; }
    }
    if (ok) bag_to_red_cls[b].push_back(cl_p);
    // if not ok, then red would contract two components, so ignore
  }

  uint32_t cls_added = 0;
  for(uint32_t i = 0; i < bags.size(); i++) {
    const auto& bag = bags[i];
    CDat sub_c;
    sub_c.conf = cdat->conf;
    sub_c.fg = cdat->fg->dup();
    assert(std::is_sorted(bag.begin(), bag.end()));

    // Assign var numbers: indeps first (1..k), then opt_indeps (k+1..k+m), then rest
    vector<int> var_map(cdat->nvars+1, -1);
    sub_c.nvars = bag.size();
    uint32_t next_num = 1;
    for(const auto& v: bag) if (cdat->indeps.count(v)) var_map[v] = next_num++;
    for(const auto& v: bag) if (!cdat->indeps.count(v) && cdat->opt_indeps.count(v)) var_map[v] = next_num++;
    for(const auto& v: bag) if (!cdat->indeps.count(v) && !cdat->opt_indeps.count(v)) var_map[v] = next_num++;
    for(const auto& v: bag) {
      if (cdat->indeps.count(v)) sub_c.indeps.insert(var_map[v]);
      if (cdat->opt_indeps.count(v)) sub_c.opt_indeps.insert(var_map[v]);
      if (auto it = cdat->lit_weights.find(GanakInt::Lit(v, false)); it != cdat->lit_weights.end())
        sub_c.lit_weights[GanakInt::Lit(var_map[v], false)] = it->second->dup();
      if (auto it = cdat->lit_weights.find(GanakInt::Lit(v, true)); it != cdat->lit_weights.end())
        sub_c.lit_weights[GanakInt::Lit(var_map[v], true)] = it->second->dup();
    }

    // remap clauses
    auto remap_clause = [&](const vector<GanakInt::Lit>& cl) {
      vector<GanakInt::Lit> new_cl;
      new_cl.reserve(cl.size());
      for(const auto& l: cl) {
        assert(var_map[l.var()] != -1);
        assert(var_map[l.var()] < (int)sub_c.nvars +1);
        new_cl.push_back(GanakInt::Lit(var_map[l.var()], l.sign()));
      }
      return new_cl;
    };
    for(const auto& cl: bag_to_irred_cls[i]) sub_c.irred_cls.emplace_back(remap_clause(cl));
    for(const auto& cl: bag_to_red_cls[i]) sub_c.red_cls.push_back({remap_clause(cl.first), cl.second});
    cls_added += sub_c.irred_cls.size();

    // set up counter
    if (cdat->conf.verb >= 2)
      cout << "c o Counting component with " << sub_c.nvars << " vars, "
        << sub_c.irred_cls.size() << " irredundant clauses, "
        << sub_c.red_cls.size() << " redundant clauses, "
        << sub_c.indeps.size() << " independent vars, "
        << sub_c.opt_indeps.size() << " optional independent vars" << endl;

    if (sub_c.indeps.size() == 0) sub_c.conf.verb = 0;
    if (sub_c.indeps.size() == 0 && sub_c.irred_cls.size() == 0) continue;
    if (sub_c.indeps.size() == 0 && sub_c.irred_cls.size() < 10) {
      assert(!sub_c.irred_cls.empty());
      bool all_same = true;
      auto one = sub_c.irred_cls[0];
      assert(!one.empty());
      for(size_t i2 = 1; i2 < sub_c.irred_cls.size(); i2++) {
        if (sub_c.irred_cls[i2] != one) { all_same = false; break; }
      }
      if (all_same) continue;
    }

    // run counter
    OuterCounter out_cnt(sub_c.conf, sub_c.fg->dup());
    out_cnt.new_vars(sub_c.nvars);
    out_cnt.set_indep_support(sub_c.indeps);
    out_cnt.set_optional_indep_support(sub_c.opt_indeps);
    for(const auto& w: sub_c.lit_weights) out_cnt.set_lit_weight(w.first, w.second);
    for(const auto& cl: sub_c.irred_cls) out_cnt.add_irred_cl(cl);
    for(const auto& p: sub_c.red_cls) out_cnt.add_red_cl(p.first, p.second);
    if (sub_c.conf.verb && cdat->print_indep_distrib) out_cnt.print_indep_distrib();
    auto ret = out_cnt.count(bits_jobs, num_threads);
    *cnt *= *ret;
    if (sub_c.conf.verb) cout << "c o intermediate count: " << *ret << endl;
    cdat->is_approximate |= out_cnt.get_is_approximate();
    cdat->max_num_cache_lookups = std::max(cdat->max_num_cache_lookups, out_cnt.get_num_cache_lookups());
    cdat->max_cache_elems = std::max(cdat->max_cache_elems, out_cnt.get_max_cache_elems());
  }
  assert(cls_added == cdat->irred_cls.size());
  return cnt;
}

/* DLL_PUBLIC void Ganak::set_generators(const std::vector<std::map<GanakInt::Lit, GanakInt::Lit>>& _gens) { */
/*   CDat* c = (CDat*)cdat; */
/*   cdat->set_generators(_gens); */
/* } */
DLL_PUBLIC void Ganak::set_indep_support(const std::set<uint32_t>& indeps) {
  for(const auto& v: indeps) {
    if (v > cdat->nvars) {
      cerr << "ERROR: setting independent support variable " << v
           << " larger than number of variables: " << cdat->nvars << endl;
      assert(false);
      exit(EXIT_FAILURE);
    }
  }
  cdat->indeps = indeps;
}
DLL_PUBLIC bool Ganak::add_red_cl(const std::vector<GanakInt::Lit>& lits, int lbd) {
  cdat->red_cls.push_back({lits, (uint32_t)lbd});
  return true;
}
DLL_PUBLIC bool Ganak::get_is_approximate() const {
  return cdat->is_approximate;
}
DLL_PUBLIC bool Ganak::add_irred_cl(const std::vector<GanakInt::Lit>& lits) {
  cdat->irred_cls.push_back(lits);
  return true;
}
DLL_PUBLIC void Ganak::set_optional_indep_support(const std::set<uint32_t>& indeps) {
  for(const auto& v: indeps) {
    if (v > cdat->nvars) {
      cerr << "ERROR: setting optional independent support variable " << v
           << " larger than number of variables: " << cdat->nvars << endl;
      assert(false);
      exit(EXIT_FAILURE);
    }
  }
  cdat->opt_indeps = indeps;
}
DLL_PUBLIC void Ganak::set_lit_weight(const GanakInt::Lit l, const FF& w) {
  if (cdat->opt_indeps.count(l.var()) == 0 && cdat->indeps.count(l.var()) == 0) {
    cerr << "ERROR: setting weight for literal " << (l.sign() ? "" : "-") << l.var()
         << " which is not in independent or optional independent support" << endl;
    assert(false);
    exit(EXIT_FAILURE);
  }
  cdat->lit_weights[l] = w->dup();
}
DLL_PUBLIC void Ganak::new_vars(const uint32_t n) {
  assert(cdat->nvars == 0);
  cdat->nvars = n;
}
DLL_PUBLIC void Ganak::print_indep_distrib() const {
  cdat->print_indep_distrib = true;
}

DLL_PUBLIC uint64_t Ganak::get_num_cache_lookups() const {
  return cdat->max_num_cache_lookups;
}

DLL_PUBLIC uint64_t Ganak::get_max_cache_elems() const {
  return cdat->max_cache_elems;
}

vector<vector<uint32_t>> find_disconnected(const CDat& dat) {
  vector<int> var_to_bag(dat.nvars+1, -1);
  map<int, vector<int>> bag_to_vars;
  int bag_id = 0;
  set<int> bags;

  auto move_to_bag = [&](int bid, int v) {
    if (var_to_bag[v] == bid) return;
    if (var_to_bag[v] == -1) {
      var_to_bag[v] = bid;
      bag_to_vars[bid].push_back(v);
      return;
    }
    assert (var_to_bag[v] != bid);
    int old_bid = var_to_bag[v];
    for(const auto& vv: bag_to_vars[old_bid]) var_to_bag[vv] = bid;
    bag_to_vars[bid].insert(bag_to_vars[bid].end(), bag_to_vars[old_bid].begin(), bag_to_vars[old_bid].end());
    bag_to_vars.erase(old_bid);
    bags.erase(old_bid);
  };

  for(const auto& cl: dat.irred_cls) {
    assert(!cl.empty() && "We filtered out empty clauses above");
    set<int> vars_in_cl;
    for(const auto& l: cl) vars_in_cl.insert(l.var());

    bool found = false;
    for(const auto& v: vars_in_cl) {
      if (var_to_bag[v] != -1) {
        found = true;
        int bid = var_to_bag[v];
        for(const auto& vv: vars_in_cl) move_to_bag(bid, vv);
        break;
      }
    }
    if (!found) {
      bag_to_vars[bag_id] = {};
      for(const auto& v: vars_in_cl) {
        var_to_bag[v] = bag_id;
        bag_to_vars[bag_id].push_back(v);
      }
      bags.insert(bag_id);
      bag_id++;
    }
  }
  for(uint32_t i = 1; i <= dat.nvars; i++) {
    if (var_to_bag[i] == -1) {
      bag_to_vars[bag_id] = {(int)i};
      var_to_bag[i] = bag_id;
      bags.insert(bag_id);
      bag_id++;
    }
  }

  vector<vector<uint32_t>> res;
  res.reserve(bags.size());
  for(const auto& b: bags) {
    /* cout << "c Found bag " << b << " with vars: "; */
    /* for(const auto& v: bag_to_vars[b]) cout << v << " "; */
    /* cout << endl; */
    const auto& bvars = bag_to_vars[b];
    res.emplace_back(bvars.begin(), bvars.end());
  }

  // Check
  uint32_t total_vars = 0;
  for(const auto& b: bags) total_vars += bag_to_vars[b].size();
  /* cout << "c Found " << bags.size() << " bags with total vars: " << total_vars << endl; */
  /* cout << "c Total vars in formula: " << dat.nvars << endl; */
  vector<int> count_vars(dat.nvars+1, 0);
  for(const auto& b: bags) {
    assert(!bag_to_vars[b].empty());
    for(const auto& v: bag_to_vars[b]) count_vars[v]++;
  }
  for(uint32_t i = 1; i <= dat.nvars; i++) {
    if (count_vars[i] != 1) {
      cout << "c ERROR: var " << i << " occurs in " << count_vars[i] << " bags" << endl;
      cout << "var_to_bag[" << i << "] = " << var_to_bag[i] << endl;
      assert(false);
    }
  }
  assert(total_vars == dat.nvars);

  if (dat.conf.verb >= 1) cout << "c o Found " << bags.size() << " disconnected component(s)" << endl;
  for(auto& b: res) std::sort(b.begin(), b.end());
  return res;
}

/* auto run_breakid(const ArjunNS::SimplifiedCNF& cnf) { */
/*   double my_time = cpu_time(); */
/*   vector<map<Lit, Lit>> generators; */
/*   BID::BreakID breakid; */
/*   /1* breakid.set_useMatrixDetection(conf.useMatrixDetection); *1/ */
/*   /1* breakid.set_useFullTranslation(conf.useFullTranslation); *1/ */
/*   breakid.set_verbosity(0); */
/*   breakid.start_dynamic_cnf(cnf.nVars()); */
/*   for(const auto& cl: cnf.clauses) { */
/*     breakid.add_clause((BID::BLit*)cl.data(), cl.size()); */
/*   } */
/*   breakid.set_steps_lim(4000); */
/*   breakid.end_dynamic_cnf(); */
/*   verb_print(1, "[breakid] Num generators: " << breakid.get_num_generators()); */
/*   breakid.detect_subgroups(); */
/*   if (conf.verb >= 1) breakid.print_generators(std::cout, "c o "); */
/*   vector<unordered_map<BID::BLit, BID::BLit> > orig_gen; */
/*   breakid.get_perms(&orig_gen); */
/*   for(const auto& m: orig_gen) { */
/*     map<Lit, Lit> gen; */
/*     for(const auto& gp: m) { */
/*       gen[Lit(gp.first.var()+1, gp.first.sign())] = Lit(gp.second.var()+1, gp.second.sign()); */
/*       if (conf.verb >= 2) cout << "c o " << gp.first << " -> " << gp.second << endl; */
/*     } */
/*     generators.push_back(gen); */
/*   } */
/*   verb_print(1, "[breakid] T: " << (cpu_time()-my_time)); */
/*   return generators; */
/* } */

