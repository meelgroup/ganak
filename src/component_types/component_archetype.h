/*
 * component_archetype.h
 *
 *  Created on: Feb 9, 2013
 *      Author: mthurley
 *  Changed in 2020 by VincentDerk, smsharma1 and tvanbr
 */

#ifndef COMPONENT_ARCHETYPE_H_
#define COMPONENT_ARCHETYPE_H_

#include "../primitive_types.h"
#include "../solver_config.h"
#include "../containers.h"
#include "component.h"
#include "cacheable_component.h"

#include "cacheable_component.h"

#include "ClHashComponent.h"
#include "sparsegraph_component.h"

// #define VERB

#include <algorithm>
#include <cstring>
#include <memory>
#include <unordered_map>


#include <iostream>
// State values for variables found during component
// analysis (CA)
typedef unsigned char CA_SearchState;
#define CA_NIL 0
#define CA_VAR_IN_SUP_COMP_UNSEEN 1
#define CA_VAR_SEEN 2
#define CA_VAR_IN_OTHER_COMP 4

#define CA_VAR_MASK 7

#define CA_CL_IN_SUP_COMP_UNSEEN 8
#define CA_CL_SEEN 16
#define CA_CL_IN_OTHER_COMP 32
#define CA_CL_ALL_LITS_ACTIVE 64

#define CA_CL_MASK 120

class StackLevel;


class ComponentArchetype
{
public:
    ComponentArchetype(SolverConfiguration &config) : config_(config)
    {
    }

    ComponentArchetype(StackLevel &stack_level, Component &super_comp, SolverConfiguration &config) : p_super_comp_(&super_comp), p_stack_level_(&stack_level), config_(config)
    {
    }

    ~ComponentArchetype() {
        delete[] seen_;
        delete[] idMap_;
    }

    void reInitialize(StackLevel &stack_level, Component &super_comp)
    {
        p_super_comp_ = &super_comp;
        p_stack_level_ = &stack_level;
        clearArrays();
        current_comp_for_caching_.reserveSpace(super_comp.num_variables(), super_comp.numLongClauses());
    }

    Component &super_comp()
    {
        return *p_super_comp_;
    }

    StackLevel &stack_level()
    {
        return *p_stack_level_;
    }

    void setVar_in_sup_comp_unseen(VariableIndex v)
    {
        seen_[v] = CA_VAR_IN_SUP_COMP_UNSEEN | (seen_[v] & CA_CL_MASK);
    }

    void setClause_in_sup_comp_unseen(ClauseIndex cl)
    {
        seen_[cl] = CA_CL_IN_SUP_COMP_UNSEEN | (seen_[cl] & CA_VAR_MASK);
    }

    void setVar_nil(VariableIndex v)
    {
        seen_[v] &= CA_CL_MASK;
    }

    void setClause_nil(ClauseIndex cl)
    {
        seen_[cl] &= CA_VAR_MASK;
    }

    void setVar_seen(VariableIndex v)
    {
        seen_[v] = CA_VAR_SEEN | (seen_[v] & CA_CL_MASK);
    }

    void setClause_seen(ClauseIndex cl)
    {
        setClause_nil(cl);
        seen_[cl] = CA_CL_SEEN | (seen_[cl] & CA_VAR_MASK);
    }

    void setClause_seen(ClauseIndex cl, bool all_lits_act)
    {
        setClause_nil(cl);
        seen_[cl] = CA_CL_SEEN | (all_lits_act ? CA_CL_ALL_LITS_ACTIVE : 0) | (seen_[cl] & CA_VAR_MASK);
    }

    void setVar_in_other_comp(VariableIndex v)
    {
        seen_[v] = CA_VAR_IN_OTHER_COMP | (seen_[v] & CA_CL_MASK);
    }

    void setClause_in_other_comp(ClauseIndex cl)
    {
        seen_[cl] = CA_CL_IN_OTHER_COMP | (seen_[cl] & CA_VAR_MASK);
    }

    bool var_seen(VariableIndex v)
    {
        return seen_[v] & CA_VAR_SEEN;
    }

    bool clause_seen(ClauseIndex cl)
    {
        return seen_[cl] & CA_CL_SEEN;
    }

    bool clause_all_lits_active(ClauseIndex cl)
    {
        return seen_[cl] & CA_CL_ALL_LITS_ACTIVE;
    }

    void setClause_all_lits_active(ClauseIndex cl)
    {
        seen_[cl] |= CA_CL_ALL_LITS_ACTIVE;
    }

    bool var_nil(VariableIndex v)
    {
        return (seen_[v] & CA_VAR_MASK) == 0;
    }

    bool clause_nil(ClauseIndex cl)
    {
        return (seen_[cl] & CA_CL_MASK) == 0;
    }

    bool var_unseen_in_sup_comp(VariableIndex v)
    {
        return seen_[v] & CA_VAR_IN_SUP_COMP_UNSEEN;
    }

    bool clause_unseen_in_sup_comp(ClauseIndex cl)
    {
        return seen_[cl] & CA_CL_IN_SUP_COMP_UNSEEN;
    }

    bool var_seen_in_peer_comp(VariableIndex v)
    {
        return seen_[v] & CA_VAR_IN_OTHER_COMP;
    }

    bool clause_seen_in_peer_comp(ClauseIndex cl)
    {
        return seen_[cl] & CA_CL_IN_OTHER_COMP;
    }

    static void initArrays(unsigned max_variable_id, unsigned max_clause_id)
    {
        unsigned seen_size = std::max(max_variable_id, max_clause_id) + 1;
        seen_ = new CA_SearchState[seen_size];
        seen_byte_size_ = sizeof(CA_SearchState) * (seen_size);
        clearArrays();
        idMap_ = new unsigned[(max_variable_id + 1) << 1];
    }

    static void clearArrays()
    {
        memset(seen_, CA_NIL, seen_byte_size_);
    }

    SparseGraphCacheComponent *constructIsoComp(Component *comp, const shared_ptr<std::vector<LiteralID>> &lit_pool,
                                                std::vector<ClauseOfs> &clauseOffsets,
                                                LiteralIndexedVector<Literal> &literals)
    {
        unsigned num_variables = comp->num_variables();
        unsigned num_clauses = comp->numLongClauses();

        int n = num_variables * 2 + num_clauses;
        auto *edges = new std::vector<unsigned>[n]; // For each vertex, all its connections.
        int num_edges = 0;

        unsigned nodeIndex = 0;
        /* Create lit - lit edges */
        for (auto v_it = comp->varsBegin(); *v_it != varsSENTINEL; v_it++)
        {
            auto lit_f = LiteralID(*v_it, false);
            edges[nodeIndex] = std::vector<unsigned>();
            edges[nodeIndex].reserve(num_clauses + 1 + literals[lit_f].binary_links_orig_size);
            auto lit_t = LiteralID(*v_it, true);
            edges[nodeIndex + 1] = std::vector<unsigned>();
            edges[nodeIndex + 1].reserve(num_clauses + literals[lit_t].binary_links_orig_size);

            edges[nodeIndex].push_back(nodeIndex + 1); //connected to negation
            edges[nodeIndex + 1].push_back(nodeIndex); //connected to negation
            num_edges++;
            //std::cout << "connecting " << nodeIndex << "(" << (unsigned) LiteralID(*v_it, false) << ") and " << (nodeIndex + 1) << "\n";
            idMap_[(unsigned)lit_f] = nodeIndex++;
            idMap_[(unsigned)lit_t] = nodeIndex++;
        }

        /* Create cl-lit edges */
        for (auto it_cl = comp->clsBegin(); *it_cl != clsSENTINEL; it_cl++)
        {
            edges[nodeIndex] = std::vector<unsigned>();
            edges[nodeIndex].reserve(num_variables);
            ClauseOfs offset = clauseOffsets[*it_cl];
            for (auto it_lit_cl = (*lit_pool).begin() + static_cast<unsigned>(offset);
                 *it_lit_cl != SENTINEL_LIT; it_lit_cl++)
            {
                if (var_seen(it_lit_cl->var()))
                {
                    unsigned litIndex = idMap_[(unsigned)*it_lit_cl];
                    edges[nodeIndex].push_back(litIndex);
                    edges[litIndex].push_back(nodeIndex);
                    num_edges++;
                    //std::cout << "connecting clause " << nodeIndex << "and " << litIndex << "\n";
                }
            }
            nodeIndex++;
        }

        /* Create cl-lit binary edges */
        for (auto v_it = comp->varsBegin(); *v_it != varsSENTINEL; v_it++)
        {
            LiteralID lit = LiteralID(*v_it, false);
            unsigned currVarNodeIndex = idMap_[(unsigned)lit];
            unsigned b_index = 0;
            for (auto b_it = literals[lit].binary_links_.begin(); b_index < literals[lit].binary_links_orig_size && *b_it != SENTINEL_LIT; b_it++, b_index++)
            {
                if (var_seen(b_it->var()))
                {
                    unsigned otherVarNodeIndex = idMap_[(unsigned)*b_it];
                    edges[currVarNodeIndex].push_back(otherVarNodeIndex);
                    edges[otherVarNodeIndex].push_back(currVarNodeIndex);
                    //std::cout << "connecting binary " << currVarNodeIndex << "and " << otherVarNodeIndex << "\n";
                    num_edges++;
                }
            }

            lit = LiteralID(*v_it, true);
            currVarNodeIndex = idMap_[(unsigned)lit];
            b_index = 0;
            for (auto b_it = literals[lit].binary_links_.begin(); b_index < literals[lit].binary_links_orig_size && *b_it != SENTINEL_LIT; b_it++, b_index++)
            {
                if (var_seen(b_it->var()))
                {
                    unsigned otherVarNodeIndex = idMap_[(unsigned)*b_it];
                    edges[currVarNodeIndex].push_back(otherVarNodeIndex);
                    edges[otherVarNodeIndex].push_back(currVarNodeIndex);
                    //std::cout << "connecting binary " << currVarNodeIndex << "and " << otherVarNodeIndex << "\n";
                    num_edges++;
                }
            }
        }
        num_edges *= 2;

        //nv = #vertices = num_clauses + num_variables * 2;
        //d = array, for each vertex the degree of edges.
        //v = array, for each vertex the start index in e
        //e = array of endpoint of edges such that e[v[i]], e[v[i]+1], ..., e[v[i]+d-1] are connected to vertex i
        DYNALLSTAT(int, lab, lab_sz);
        DYNALLSTAT(int, ptn, ptn_sz);
        DYNALLSTAT(int, orbits, orbits_sz);
        static DEFAULTOPTIONS_SPARSEGRAPH(options);
        options.getcanon = TRUE;
        statsblk stats;
        SG_DECL(sg);

        int m = SETWORDSNEEDED(n);
        nauty_check(WORDSIZE, m, n, NAUTYVERSIONID);
        DYNALLOC1(int, lab, lab_sz, n, "malloc");
        DYNALLOC1(int, ptn, ptn_sz, n, "malloc");
        DYNALLOC1(int, orbits, orbits_sz, n, "malloc");
        SG_ALLOC(sg, n, num_edges, "malloc");

        /* Create colour partitions */
        for (int i = 0; i < n; i++)
        {
            lab[i] = i;
            ptn[i] = 1;
        }
        ptn[num_variables * 2 - 1] = 0;

        /* Create graph */
        //sg.d[i] = degree of node i
        //sg.v[i] = index in sg.e from which all neighbors of i are listed (sg.d[i] neighbors).
        //sg.e = concatenation of all neighbors of node 0, node 1, node 2, ...
        sg.nv = n;
        sg.nde = num_edges;

        //Fill in sg.d, sg.v and sg.e.
        // for index 0 (separate because sg.v[-1] does not exist)
        sg.d[0] = edges[0].size();
        sg.v[0] = 0;
        std::memcpy(sg.e, edges[0].data(), sizeof(unsigned) * edges[0].size());
        // for index 1..n
        for (int i = 1; i < n; i++)
        {
            sg.d[i] = edges[i].size();
            sg.v[i] = sg.v[i - 1] + edges[i - 1].size();
            std::memcpy(&(sg.e[sg.v[i]]), edges[i].data(), sizeof(unsigned) * edges[i].size());
        }

        //std::cout << "Start nauty with num_var: " << num_variables << ", num_clauses " << num_clauses << " and num_edges " << num_edges << "\n";
        auto *cg = (sparsegraph *)malloc(sizeof(sparsegraph));
        SG_INIT(*cg);
        sparsenauty(&sg, lab, ptn, orbits, &options, &stats, cg);
        long cg_hashkey = hashgraph_sg(cg, 0);
        SG_FREE(sg);
        delete[] edges;

        auto new_comp = new SparseGraphCacheComponent(cg, cg_hashkey, num_variables);
        if (config_.use_icsvsads)
        {
            for (auto it = comp->varsBegin(); *it != varsSENTINEL; it++)
            {
                new_comp->variables_.insert(*it);
            }
        }
        return new_comp;
    }

    /**
     * Construct a component holding the CLhash of the STD encoding
     * Beware, this can only be called with a component that was only just obtained from makeComponentFromState().
     * After this, setComponentVarsSeen() should be called:
     * comp = makeComponentFromState()
     * cg = constructSTDClHashComp(...)
     * setComponentVarsSeen(comp)
     * @param comp The component for which to construct the Cl Hash Comp. The variables in comp must be exactly the
     * variables regarded as part of the current state. (cf. setComponentVarsSeen() and makeComponentFromState)!
     * @param lit_pool
     * @param clauseOffsets
     * @return
     */
    ClHashComponent *constructSTDClHashComp(vector<void *> &random, Component *comp, shared_ptr<std::vector<LiteralID>> lit_pool,
                                            std::vector<ClauseOfs> &clauseOffsets)
    {
        std::vector<unsigned> data_vector; //TODO: Figure out a decent space reserve?

        //For each clause of comp, store the literals that are active (separate clause with varsSENTINEL).
        for (auto it_cl = comp->clsBegin(); *it_cl != clsSENTINEL; it_cl++)
        {
            ClauseOfs offset = clauseOffsets[*it_cl - 1];
            for (auto it_lit_cl = (*lit_pool).begin() + static_cast<unsigned>(offset);
                 *it_lit_cl != SENTINEL_LIT; it_lit_cl++)
            {
                if (var_seen(it_lit_cl->var()))
                    data_vector.push_back((unsigned)*it_lit_cl);
            }
            data_vector.push_back(varsSENTINEL);
        }

        //Create clhash
        auto *clhashkey = new uint64_t[random.capacity()];
        for (unsigned int i = 0; i < random.capacity(); i++)
        {
            clhasher h(random[i]);
            clhashkey[i] = h(data_vector);
        }

        return new ClHashComponent(clhashkey, random.capacity(), comp->num_variables());
    }

    /**
     * Set variables of the given component as seen (set as part of other component).
     * @param comp The component with variables to set as seen.
     */
    void setComponentVarsSeen(Component *comp)
    {
        for (auto v_it = comp->varsBegin(); *v_it != varsSENTINEL; v_it++)
            setVar_in_other_comp(*v_it);
    }

    /**
     * Make a component from the current state and set the variables of this component as seen
     * (='part of other component'). The variables are not set as seen, execute setComponentVarsSeen(comp) after this.
     * This allows us to use the SEEN bitset even after calling this, required for constructSparseCG.
     * @param stack_size The amount of variables to make space for
     * @return The component associated with the current state
     */
    Component *makeComponentFromState(unsigned stack_size)
    {
        auto *p_new_comp = new Component();
        p_new_comp->reserveSpace(stack_size, super_comp().numLongClauses());
        current_comp_for_caching_.clear();

        for (auto v_it = super_comp().varsBegin(); *v_it != varsSENTINEL; v_it++)
            if (var_seen(*v_it))
            { //we have to put a var into our component
                p_new_comp->addVar(*v_it);
                current_comp_for_caching_.addVar(*v_it);
                //setVar_in_other_comp(*v_it);
            }
        p_new_comp->closeVariableData();
        current_comp_for_caching_.closeVariableData();

        for (auto it_cl = super_comp().clsBegin(); *it_cl != clsSENTINEL; it_cl++)
            if (clause_seen(*it_cl))
            {
                p_new_comp->addCl(*it_cl);
                if (!clause_all_lits_active(*it_cl))
                    current_comp_for_caching_.addCl(*it_cl);
                setClause_in_other_comp(*it_cl);
            }
        p_new_comp->closeClauseData();
        current_comp_for_caching_.closeClauseData();
        return p_new_comp;
    }


    Component current_comp_for_caching_;

private:
    Component *p_super_comp_;
    StackLevel *p_stack_level_;
    SolverConfiguration &config_;

    static CA_SearchState *seen_;
    static unsigned seen_byte_size_;
    static unsigned *idMap_; //Used in Canonical Graph creation, mapping literalIDs to graph index
};

#endif /* COMPONENT_ARCHETYPE_H_ */
