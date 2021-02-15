/*
 * solver.h
 *
 *  Created on: Aug 23, 2012
 *      Author: marc
 */

#ifndef SOLVER_H_
#define SOLVER_H_


#include "statistics.h"
#include "instance.h"
#include "component_management.h"



#include "solver_config.h"

#include <sys/time.h>

enum retStateT {
	EXIT, RESOLVED, PROCESS_COMPONENT, BACKTRACK, RESTART
};



class StopWatch {
public:

  StopWatch();

  bool timeBoundBroken() {
    timeval actual_time;
    gettimeofday(&actual_time, NULL);
    return actual_time.tv_sec - start_time_.tv_sec > time_bound_;
  }

  bool start() {
    bool ret = gettimeofday(&last_interval_start_, NULL);
    start_time_ = stop_time_ = last_interval_start_;
    return !ret;
  }

  bool stop() {
    return gettimeofday(&stop_time_, NULL) == 0;
  }

  double getElapsedSeconds() {
    timeval r = getElapsedTime();
    return r.tv_sec + (double) r.tv_usec / 1000000;
  }

  bool interval_tick() {
    timeval actual_time;
    gettimeofday(&actual_time, NULL);
    if (actual_time.tv_sec - last_interval_start_.tv_sec
        > interval_length_.tv_sec) {
      gettimeofday(&last_interval_start_, NULL);
      return true;
    }
    return false;
  }

  void setTimeBound(long int seconds) {
    time_bound_ = seconds;
  }
  long int getTimeBound();

private:
  timeval start_time_;
  timeval stop_time_;

  long int time_bound_;

  timeval interval_length_;
  timeval last_interval_start_;
	

  // if we have started and then stopped the watch, this returns
  // the elapsed time
  // otherwise, time elapsed from start_time_ till now is returned
  timeval getElapsedTime();
};

class Solver: public Instance {
public:
	Solver() {
	  stopwatch_.setTimeBound(config_.time_bound_seconds);
		literal_pool_ = make_shared<vector<LiteralID>> ();
	}

	void solve(const string & file_name);

	SolverConfiguration &config() {
		return config_;
	}

	DataAndStatistics &statistics() {
	        return statistics_;
	}
	void setTimeBound(long int i) {
		stopwatch_.setTimeBound(i);
	}

private:
	SolverConfiguration config_;
	bool isindependent = true;

	DecisionStack stack_; // decision stack
	vector<LiteralID> literal_stack_;

	StopWatch stopwatch_;

	ComponentManager comp_manager_ =
	  ComponentManager(
	    config_,
	    statistics_,
		literal_values_,
		literals_);

	// the last time conflict clauses have been deleted
	unsigned long last_ccl_deletion_time_ = 0;
	// the last time the conflict clause storage has been compacted
	unsigned long last_ccl_cleanup_time_ = 0;

	bool simplePreProcess();
	bool prepFailedLiteralTest();
	// we assert that the formula is consistent
	// and has not been found UNSAT yet
	// hard wires all assertions in the literal stack into the formula
	// removes all set variables and essentially reinitiallizes all
	// further data
	void HardWireAndCompact();

	SOLVER_StateT countSAT();

	void decideLiteral();
	bool bcp();


	 void decayActivitiesOf(Component & comp) {
	   for (auto it = comp.varsBegin(); *it != varsSENTINEL; it++) {
	          literal(LiteralID(*it,true)).activity_score_ *=0.5;
	          literal(LiteralID(*it,false)).activity_score_ *=0.5;
	       }
	}
	///  this method performs Failed literal tests online
	bool implicitBCP();

	// this is the actual BCP algorithm
	// starts propagating all literal in literal_stack_
	// beginingg at offset start_at_stack_ofs
	bool BCP(unsigned start_at_stack_ofs);

	retStateT backtrack();

	// if on the current decision level
	// a second branch can be visited, RESOLVED is returned
	// otherwise returns BACKTRACK
	retStateT resolveConflict();

	/////////////////////////////////////////////
	//  BEGIN small helper functions
	/////////////////////////////////////////////

	float scoreOf(VariableIndex v) {
		float score = comp_manager_.scoreOf(v);
		score += 10.0 * literal(LiteralID(v, true)).activity_score_;
		score += 10.0 * literal(LiteralID(v, false)).activity_score_;
//		score += (10*stack_.get_decision_level()) * literal(LiteralID(v, true)).activity_score_;
//		score += (10*stack_.get_decision_level()) * literal(LiteralID(v, false)).activity_score_;

		return score;
	}

	bool setLiteralIfFree(LiteralID lit,
			Antecedent ant = Antecedent(NOT_A_CLAUSE)) {
		if (literal_values_[lit] != X_TRI)
			return false;
		var(lit).decision_level = stack_.get_decision_level();
		var(lit).ante = ant;
		var(lit).polarity = lit.sign();
		var(lit).set = true;
		literal_stack_.push_back(lit);
		if (ant.isAClause() && ant.asCl() != NOT_A_CLAUSE)
			getHeaderOf(ant.asCl()).increaseScore();
		literal_values_[lit] = T_TRI;
		literal_values_[lit.neg()] = F_TRI;
		return true;
	}

	void printOnlineStats();

	void print(vector<LiteralID> &vec);
	void print(vector<unsigned> &vec);


	void setConflictState(LiteralID litA, LiteralID litB) {
		violated_clause.clear();
		violated_clause.push_back(litA);
		violated_clause.push_back(litB);
	}
	void setConflictState(ClauseOfs cl_ofs) {
		getHeaderOf(cl_ofs).increaseScore();
		violated_clause.clear();
		for (auto it = beginOf(cl_ofs); *it != SENTINEL_LIT; it++)
			violated_clause.push_back(*it);
	}

	vector<LiteralID>::const_iterator TOSLiteralsBegin() {
		return literal_stack_.begin() + stack_.top().literal_stack_ofs();
	}

	void initStack(unsigned int resSize) {
		stack_.clear();
		stack_.reserve(resSize);
		literal_stack_.clear();
		literal_stack_.reserve(resSize);
		// initialize the stack to contain at least level zero
		stack_.push_back(StackLevel(1, 0, 2));
		stack_.back().changeBranch();
	}

	const LiteralID &TOS_decLit() {
		assert(stack_.top().literal_stack_ofs() < literal_stack_.size());
		return literal_stack_[stack_.top().literal_stack_ofs()];
	}

	void reactivateTOS() {
		for (auto it = TOSLiteralsBegin(); it != literal_stack_.end(); it++)
			unSet(*it);
		comp_manager_.cleanRemainingComponentsOf(stack_.top());
		literal_stack_.resize(stack_.top().literal_stack_ofs());
		stack_.top().resetRemainingComps();
	}

	bool fail_test(LiteralID lit) {
		unsigned sz = literal_stack_.size();
		// we increase the decLev artificially
		// s.t. after the tentative BCP call, we can learn a conflict clause
		// relative to the assignment of *jt
		stack_.startFailedLitTest();
		setLiteralIfFree(lit);

		assert(!hasAntecedent(lit));

		bool bSucceeded = BCP(sz);
		if (!bSucceeded)
			recordAllUIPCauses();

		stack_.stopFailedLitTest();

		while (literal_stack_.size() > sz) {
			unSet(literal_stack_.back());
			literal_stack_.pop_back();
		}
		return bSucceeded;
	}
	/////////////////////////////////////////////
	//  BEGIN conflict analysis
	/////////////////////////////////////////////

	// if the state name is CONFLICT,
	// then violated_clause contains the clause determining the conflict;
	vector<LiteralID> violated_clause;
	// this is an array of all the clauses found
	// during the most recent conflict analysis
	// it might contain more than 2 clauses
	// but always will have:
	//      uip_clauses_.front() the 1UIP clause found
	//      uip_clauses_.back() the lastUIP clause found
	//  possible clauses in between will be other UIP clauses
	vector<vector<LiteralID> > uip_clauses_;

	// the assertion level of uip_clauses_.back()
	// or (if the decision variable did not have an antecedent
	// before) then assertionLevel_ == DL;
	int assertion_level_ = 0;

	// build conflict clauses from most recent conflict
	// as stored in state_.violated_clause
	// solver state must be CONFLICT to work;
	// this first method record only the last UIP clause
	// so as to create clause that asserts the current decision
	// literal
	void recordLastUIPCauses();
	void recordAllUIPCauses();

	void minimizeAndStoreUIPClause(LiteralID uipLit,
			vector<LiteralID> & tmp_clause, bool seen[]);
	void storeUIPClause(LiteralID uipLit, vector<LiteralID> & tmp_clause);
	int getAssertionLevel() const {
		return assertion_level_;
	}

	/////////////////////////////////////////////
	//  END conflict analysis
	/////////////////////////////////////////////
};

#endif /* SOLVER_H_ */
