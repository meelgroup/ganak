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
#include "counter.hpp"
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
  vector<std::pair<vector<GanakInt::Lit>, uint32_t>> red_cls;
  bool print_indep_distrib = false;
  bool is_approximate = false;
  uint64_t max_num_cache_lookups = 0;
  uint64_t max_cache_elems = 0;
};

vector<vector<uint32_t>> find_disconnected(const CDat& dat);

DLL_PUBLIC Ganak::Ganak(CounterConfiguration& conf, FG& fg) {
  cdat = new CDat;
  CDat* c = (CDat*)cdat;
  c->conf = conf;
  c->fg = fg->dup();
}

DLL_PUBLIC Ganak::~Ganak() {
  CDat* c = (CDat*)cdat;
  delete c;
  cdat = nullptr;
}

DLL_PUBLIC FF Ganak::count() {
  CDat* c = (CDat*)cdat;
  auto cnt = c->fg->one();

  // Check for empty clause
  for(const auto& cl: c->irred_cls) {
    if (cl.size() == 0) {
      cout << "c o intermediate count: " << *c->fg->zero() << endl;
      return c->fg->zero();
    }
  }

  auto bags = find_disconnected(*c);
  uint32_t cls_added = 0;
  for(auto& bag: bags) {
    vector<int> var_map(c->nvars+1, -1);
    CDat sub_c;
    sub_c.conf = c->conf;
    sub_c.fg = c->fg->dup();
    sub_c.nvars = bag.size();
    std::sort(bag.begin(), bag.end());
    for(size_t i = 0; i < bag.size(); i++) var_map[bag[i]] = i+1;
    for(const auto& v: bag) {
      if (c->indeps.count(v)) sub_c.indeps.insert(var_map[v]);
      if (c->opt_indeps.count(v)) sub_c.opt_indeps.insert(var_map[v]);
      if (c->lit_weights.count(GanakInt::Lit(v, false))) sub_c.lit_weights[GanakInt::Lit(var_map[v], false)] = c->lit_weights[GanakInt::Lit(v, false)]->dup();
      if (c->lit_weights.count(GanakInt::Lit(v, true))) sub_c.lit_weights[GanakInt::Lit(var_map[v], true)] = c->lit_weights[GanakInt::Lit(v, true)]->dup();
    }
    for(const auto& cl: c->irred_cls) {
      bool in_bag = true;
      vector<GanakInt::Lit> new_cl;
      for(const auto& l: cl) {
        if (var_map[l.var()] == -1) { in_bag = false; break; }
        new_cl.push_back(GanakInt::Lit(var_map[l.var()], l.sign()));
      }
      if (in_bag) sub_c.irred_cls.push_back(new_cl);
    }
    for(const auto& cl: c->red_cls) {
      bool in_bag = true;
      vector<GanakInt::Lit> new_cl;
      for(const auto& l: cl.first) {
        if (var_map[l.var()] == -1) { in_bag = false; break; }
        new_cl.push_back(GanakInt::Lit(var_map[l.var()], l.sign()));
      }
      if (in_bag) sub_c.red_cls.push_back({new_cl, cl.second});
    }
    cls_added += sub_c.irred_cls.size();

    // Now count
    OuterCounter counter(sub_c.conf, sub_c.fg->dup());
    counter.new_vars(sub_c.nvars);
    counter.set_indep_support(sub_c.indeps);
    counter.set_optional_indep_support(sub_c.opt_indeps);
    for(const auto& w: sub_c.lit_weights) counter.set_lit_weight(w.first, w.second);
    for(const auto& cl: sub_c.irred_cls) counter.add_irred_cl(cl);
    counter.end_irred_cls();
    for(const auto& p: sub_c.red_cls) counter.add_red_cl(p.first, p.second);
    if (c->print_indep_distrib) counter.print_indep_distrib();
    *cnt *= *counter.outer_count();
    cout << "c o intermediate count: " << *cnt << endl;
    c->is_approximate |= counter.get_is_approximate();
    c->max_num_cache_lookups = std::max(c->max_num_cache_lookups, counter.get_num_cache_lookups());
    c->max_cache_elems = std::max(c->max_cache_elems, counter.get_max_cache_elems());
  }
  assert(cls_added == c->irred_cls.size());
  return cnt;
}

/* DLL_PUBLIC void Ganak::set_generators(const std::vector<std::map<GanakInt::Lit, GanakInt::Lit>>& _gens) { */
/*   CDat* c = (CDat*)cdat; */
/*   c->set_generators(_gens); */
/* } */
DLL_PUBLIC void Ganak::set_indep_support(const std::set<uint32_t>& indeps) {
  CDat* c = (CDat*)cdat;
  c->indeps = indeps;
}
DLL_PUBLIC bool Ganak::add_red_cl(const std::vector<GanakInt::Lit>& lits, int lbd) {
  CDat* c = (CDat*)cdat;
  c->red_cls.push_back({lits, (uint32_t)lbd});
  return true;
}
DLL_PUBLIC bool Ganak::get_is_approximate() const {
  CDat* c = (CDat*)cdat;
  return c->is_approximate;
}
DLL_PUBLIC bool Ganak::add_irred_cl(const std::vector<GanakInt::Lit>& lits) {
  CDat* c = (CDat*)cdat;
  c->irred_cls.push_back(lits);
  return true;
}
DLL_PUBLIC void Ganak::set_optional_indep_support(const std::set<uint32_t>& indeps) {
  CDat* c = (CDat*)cdat;
  c->opt_indeps = indeps;
}
DLL_PUBLIC void Ganak::set_lit_weight(const GanakInt::Lit l, const FF& w) {
  CDat* c = (CDat*)cdat;
  c->lit_weights[l] = w->dup();
}
DLL_PUBLIC void Ganak::new_vars(const uint32_t n) {
  CDat* c = (CDat*)cdat;
  assert(c->nvars == 0);
  c->nvars = n;
}
DLL_PUBLIC void Ganak::print_indep_distrib() const {
  CDat* c = (CDat*)cdat;
  c->print_indep_distrib = true;
}

DLL_PUBLIC uint64_t Ganak::get_num_cache_lookups() const {
  CDat* c = (CDat*)cdat;
  return c->max_num_cache_lookups;
}

DLL_PUBLIC uint64_t Ganak::get_max_cache_elems() const {
  CDat* c = (CDat*)cdat;
  return c->max_cache_elems;
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
    if (cl.size() == 0) continue;
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
      bag_to_vars[bag_id] = vector<int>();
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
      bag_to_vars[bag_id] = vector<int>();
      bag_to_vars[bag_id].push_back(i);
      var_to_bag[i] = bag_id;
      bags.insert(bag_id);
      bag_id++;
    }
  }

  vector<vector<uint32_t>> res;
  for(const auto& b: bags) {
    /* cout << "c Found bag " << b << " with vars: "; */
    /* for(const auto& v: bag_to_vars[b]) cout << v << " "; */
    /* cout << endl; */
    vector<uint32_t> bag_vars;
    for(const auto& v: bag_to_vars[b]) bag_vars.push_back(v);
    res.push_back(bag_vars);
  }

  // Check
  uint32_t total_vars = 0;
  for(const auto& b: bags) total_vars += bag_to_vars[b].size();
  /* cout << "c Found " << bags.size() << " bags with total vars: " << total_vars << endl; */
  /* cout << "c Total vars in formula: " << dat.nvars << endl; */
  vector<int> count_vars(dat.nvars+1, 0);
  for(const auto& b: bags) {
    assert(bag_to_vars[b].size() > 0);
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

