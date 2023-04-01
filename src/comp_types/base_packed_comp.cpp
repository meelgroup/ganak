/*
 * base_packed_comp.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: mthurley
 */
#include "base_packed_comp.h"
#include <math.h>
#include <iostream>

BPCSizes BasePackedComponent::calcPackSize(uint32_t maxVarId, uint32_t maxClId) {

  BPCSizes sz;

  if (maxVarId == 0){
    sz.bits_per_variable = 1;
  }
  else{
    sz.bits_per_variable =  log2(maxVarId) + 1;
  }
  if(maxClId == 0){
    sz.bits_per_clause = 1;
  }
  else{
    sz.bits_per_clause   = log2(maxClId) + 1;
  }

  if (maxVarId == 0 && maxClId == 0){
    sz.bits_of_data_size = 1;
  }
  else{
    sz.bits_of_data_size = log2(maxVarId + maxClId) + 1;
  }
  assert(sz.bits_of_data_size < 32 && "Otherwise, we have an issue with nVars() that seems to convert things to uint64_t which is not allocated -- uint32_t is allocated");

  sz.variable_mask = sz.clause_mask = sz.data_size_mask = 0;
  for (uint32_t i = 0; i < sz.bits_per_variable; i++)
    sz.variable_mask = (sz.variable_mask << 1) + 1;
  for (uint32_t i = 0; i < sz.bits_per_clause; i++)
    sz.clause_mask = (sz.clause_mask << 1) + 1;
  for (uint32_t i = 0; i < sz.bits_of_data_size; i++)
    sz.data_size_mask = (sz.data_size_mask << 1) + 1;
  return sz;
}
