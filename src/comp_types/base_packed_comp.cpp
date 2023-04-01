/*
 * base_packed_comp.cpp
 *
 *  Created on: Feb 4, 2013
 *      Author: mthurley
 */
#include "base_packed_comp.h"
#include <math.h>
#include <iostream>

void BasePackedComponent::adjustPackSize(uint32_t maxVarId,
    uint32_t maxClId) {

  if (maxVarId == 0){
    _bits_per_variable = 1;
  }
  else{
    _bits_per_variable =  log2(maxVarId) + 1;
  }
  if(maxClId == 0){
    _bits_per_clause = 1;
  }
  else{
    _bits_per_clause   = log2(maxClId) + 1;
  }

  if (maxVarId == 0 && maxClId == 0){
    _bits_of_data_size = 1;
  }
  else{
    _bits_of_data_size = log2(maxVarId + maxClId) + 1;
  }
  assert(_bits_of_data_size < 32 && "Otherwise, we have an issue with nVars() that seems to convert things to uint64_t which is not allocated -- uint32_t is allocated");

  _variable_mask = _clause_mask = _data_size_mask = 0;
  for (uint32_t i = 0; i < _bits_per_variable; i++)
    _variable_mask = (_variable_mask << 1) + 1;
  for (uint32_t i = 0; i < _bits_per_clause; i++)
    _clause_mask = (_clause_mask << 1) + 1;
  for (uint32_t i = 0; i < _bits_of_data_size; i++)
    _data_size_mask = (_data_size_mask << 1) + 1;
}




