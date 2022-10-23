/*
 * component_analyzer.h
 *
 *  Created on: Feb 7, 2013
 *      Author: mthurley
 */

#ifndef COMPONENT_ANALYZER_H_
#define COMPONENT_ANALYZER_H_



#include "statistics.h"
#include "component_types/component.h"
#include "component_types/base_packed_component.h"
#include "component_types/component_archetype.h"

#include <vector>
#include <cmath>
#include <gmpxx.h>
#include "containers.h"
#include "stack.h"

struct CAClauseHeader {
  unsigned clause_id = 0;
  LiteralID lit_A;
  LiteralID lit_B;

  static unsigned overheadInLits() {
    return (sizeof(CAClauseHeader) - 2 * sizeof(LiteralID)) / sizeof(LiteralID);
  }
};

#endif /* COMPONENT_ANALYZER_H_ */
