/*
 * primitive_types.h
 *
 *  Created on: Feb 5, 2013
 *      Author: mthurley
 */

#ifndef PRIMITIVE_TYPES_H_
#define PRIMITIVE_TYPES_H_

#define varsSENTINEL  0
#define clsSENTINEL   NOT_A_CLAUSE


typedef unsigned VariableIndex;
typedef unsigned ClauseIndex;
typedef unsigned ClauseOfs;

typedef unsigned CacheEntryID;

static const ClauseIndex NOT_A_CLAUSE(0);
#define SENTINEL_CL NOT_A_CLAUSE


enum SOLVER_StateT {

  NO_STATE, SUCCESS, TIMEOUT, ABORTED
};


#ifdef DEBUG
#define toDEBUGOUT(X) cout << X;
#else
#define toDEBUGOUT(X)
#endif




#endif /* PRIMITIVE_TYPES_H_ */
