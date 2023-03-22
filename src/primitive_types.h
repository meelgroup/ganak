/*
 * primitive_types.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#pragma once

#include <cstdint>

constexpr uint32_t varsSENTINEL = 0;
constexpr uint32_t clsSENTINEL = 0;
constexpr uint32_t NOT_A_CLAUSE = 0;
constexpr uint32_t SENTINEL_CL = 0;

typedef uint32_t VariableIndex;
typedef uint32_t ClauseIndex;
typedef uint32_t ClauseOfs;
typedef uint32_t CacheEntryID;

enum SOLVER_StateT {
  NO_STATE, SUCCESS
};
