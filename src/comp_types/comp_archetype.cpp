/*
 * comp_archetype.h
 *
 *  Created on: Feb 9, 2013
 *      Author: mthurley
 */

#include "comp_archetype.hpp"

// Note that this will ensure that variables and clauses are ordered
// this means that their hash is easily comparable
template<typename T>
Comp* CompArchetype<T>::make_comp(const uint32_t comp_vars_size) {
  debug_print(COLREDBG << __PRETTY_FUNCTION__ << " start.");
  Comp *p_new_comp = new Comp();
  p_new_comp->reserve_space(comp_vars_size, super_comp().num_long_cls());

  // Fill variables in new comp
  all_vars_in_comp(super_comp(), v_it)
    if (var_visited(*v_it)) {
      p_new_comp->add_var(*v_it);
      set_var_in_peer_comp(*v_it);
    }
  p_new_comp->close_vars_data();

  // Fill (long) clause IDs in new comp
  all_cls_in_comp(super_comp(), it_cl)
    if (clause_visited(*it_cl)) {
      p_new_comp->add_cl(*it_cl);
      set_clause_in_peer_comp(*it_cl);
    }
  p_new_comp->close_cls_data();

  debug_print(COLREDBG << __PRETTY_FUNCTION__ << " finish." <<
      " New comp vars: " << p_new_comp->nVars() <<
      " long cls:" << p_new_comp->num_long_cls());
  return p_new_comp;
}

template class CompArchetype<mpz_class>;
template class CompArchetype<mpfr::mpreal>;
template class CompArchetype<mpq_class>;
