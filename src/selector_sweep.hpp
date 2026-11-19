// Recognize the complete CNF before running an exact grouped selector count.
#pragma once
#include "selector_search.hpp"

namespace GanakInt::selector {
class Encoding {
  unsigned vars;
  std::vector<Scope> positive, forbidden, targets;

public:
  explicit Encoding(unsigned n) : vars(n), targets(n+1) {}
  // 8 blocks of 32 pairs plus 20 selectors per block. Larger roots decline
  // before collecting any clauses. Every accepted clause has at most 32 literals.
  static bool fits(unsigned n, size_t clauses) { return n <= 672 && clauses <= 6800; }

  bool add(const std::vector<int>& clause) {
    if (clause.empty() || clause.size() > 32) return false;
    Scope pos, neg;
    for (int lit : clause) {
      if (!lit || lit < -int(vars) || lit > int(vars)) return false;
      (lit > 0 ? pos : neg).push_back(std::abs(lit));
    }
    if (neg.empty()) positive.push_back(std::move(pos));
    else if (pos.empty()) forbidden.push_back(std::move(neg));
    else if (pos.size() == 1 && neg.size() == 1) targets[neg[0]].push_back(pos[0]);
    else return false;
    return true;
  }

  std::optional<mpz_class> solve(RootStats* stats = nullptr, Word work_limit = 250000000) const {
    // Mixed binary clauses identify the selector variables. All other variables
    // must occur in one disjoint positive pair. Every selector must occur in
    // exactly one all-selector clause, which requires at least one true selector.
    std::vector<unsigned> pair(vars+1), side(vars+1), appearances(vars+1);
    std::vector<Scope> rows;
    unsigned pairs = 0;
    for (const auto& clause : positive) {
      unsigned selectors = 0;
      for (unsigned v : clause) selectors += !targets[v].empty();
      if (selectors == clause.size()) {
        rows.push_back(clause);
        for (unsigned v : clause) if (++appearances[v] > 1) return {};
      } else {
        if (selectors || clause.size() != 2) return {};
        unsigned a = clause[0], b = clause[1];
        if (a == b || pair[a] || pair[b]) return {};
        pair[a] = pair[b] = ++pairs;
        side[b] = 1;
      }
    }
    if (rows.empty() || rows.size() > 32) return {};
    const unsigned blocks = rows.front().size();
    if (blocks < 2 || blocks > 8 || (rows.size() > 20 && blocks > 3)) return {};
    for (unsigned v = 1; v <= vars; ++v)
      if (targets[v].empty() ? !pair[v] : appearances[v] != 1) return {};

    auto pattern = [&](const Scope& literals, Scope& scope, Point& value) {
      if (literals.size() < 2 || literals.size() > 32) return false;
      Scope ranked;
      for (unsigned v : literals) {
        if (!pair[v]) return false;
        ranked.push_back(2*pair[v] + side[v]);
      }
      std::sort(ranked.begin(), ranked.end());
      for (unsigned j = 0; j < ranked.size(); ++j) {
        unsigned p = ranked[j] / 2;
        if (j && p == scope.back()) return false;
        scope.push_back(p);
        if (ranked[j] % 2) value |= Point{1} << j;
      }
      return true;
    };
    std::map<Scope, unsigned> scopes;
    std::vector<bool> used(pairs+1);
    std::vector<unsigned> block_of(vars+1);
    std::vector<Point> point_of(vars+1);
    std::vector<Block> factors(blocks);
    for (unsigned v = 1; v <= vars; ++v) {
      if (targets[v].empty()) continue;
      Scope scope;
      if (!pattern(targets[v], scope, point_of[v])) return {};
      auto [it, fresh] = scopes.try_emplace(scope, scopes.size());
      if (fresh) {
        if (scopes.size() > blocks) return {};
        for (unsigned p : scope) {
          if (used[p]) return {};
          used[p] = true;
        }
      }
      block_of[v] = it->second;
    }
    if (scopes.size() != blocks || std::count(used.begin(), used.end(), true) != pairs) return {};
    const unsigned width = scopes.begin()->first.size();
    for (const auto& entry : scopes) if (entry.first.size() != width) return {};
    for (const auto& row : rows) {
      if (row.size() != blocks) return {};
      unsigned seen = 0;
      for (unsigned v : row) {
        unsigned b = block_of[v];
        if (seen & (1u << b)) return {};
        seen |= 1u << b;
        factors[b].required.push_back(point_of[v]);
      }
    }
    for (const auto& clause : forbidden) {
      Scope scope;
      Point value = 0;
      if (!pattern(clause, scope, value)) return {};
      auto it = scopes.find(scope);
      if (it == scopes.end()) return {};
      factors[it->second].forbidden.push_back(value);
    }
    for (auto& f : factors) {
      // Pair orientation and order are arbitrary. Align columns using selector
      // rows before comparing whole blocks and sharing their transformed tables.
      const Point anchor = f.required.front();
      std::vector<std::pair<Word, unsigned>> columns;
      for (unsigned j = 0; j < width; ++j) {
        Word signature = 0;
        for (unsigned r = 0; r < rows.size(); ++r)
          signature |= Word(((f.required[r] ^ anchor) >> j) & 1) << r;
        columns.emplace_back(signature, j);
      }
      std::sort(columns.begin(), columns.end());
      auto align = [&](Point p) {
        Point result = 0;
        for (unsigned j = 0; j < width; ++j)
          result |= (((p ^ anchor) >> columns[j].second) & 1) << j;
        return result;
      };
      for (auto& p : f.required) p = align(p);
      for (auto& p : f.forbidden) p = align(p);
      auto& fs = f.forbidden;
      std::sort(fs.begin(), fs.end());
      fs.erase(std::unique(fs.begin(), fs.end()), fs.end());
      if (fs.size() > 128) return {};
    }

    Budget budget{work_limit};
    try {
      if (rows.size() > 20) {
        WitnessSearch search(width, factors, budget);
        auto answer = search.count(rows.size());
        if (stats) *stats = {blocks, width, unsigned(rows.size()), work_limit-budget.left, true};
        return answer;
      }
      std::map<Block, unsigned> ids;
      std::vector<std::vector<Word>> tables;
      std::vector<unsigned> layout;
      for (const auto& f : factors) {
        auto [it, fresh] = ids.try_emplace(f, tables.size());
        if (fresh) tables.push_back(block_table(width, f, budget));
        layout.push_back(it->second);
      }
      // OR-convolution at the full selector mask by subset inversion. Selectors
      // are optional witnesses, so all nonempty choices in each row are counted.
      budget.spend(tables.front().size() * layout.size());
      mpz_class answer = 0, term, value;
      for (unsigned s = 0; s < tables.front().size(); ++s) {
        term = 1;
        for (unsigned b : layout) {
          const Word entry = tables[b][s];
          mpz_import(value.get_mpz_t(), 1, 1, sizeof(entry), 0, 0, &entry);
          term *= value;
        }
        if ((rows.size()-std::popcount(s)) % 2) answer -= term;
        else answer += term;
      }
      if (stats) *stats = {blocks, width, unsigned(rows.size()), work_limit-budget.left};
      return answer;
    } catch (const LimitReached&) {
      return {};
    }
  }
};
}
