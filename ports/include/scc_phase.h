#ifndef SCC_PHASE_H
#define SCC_PHASE_H

#include "graph_primitives.h"
#include "phase_conv.h"
#include <symengine/symbol.h>
#include <symengine/expression.h>
#include <vector>

struct SCCPhase {
    unsigned n_caps, n_loads, n_on_sw, n_off_sw;
    SymEngine::DenseMatrix inc_on_sw, inc_on_conv, inc_on_conv_sw, graph;
    std::vector<unsigned> sw_idxs;
    std::vector<unsigned> tree;
    SymEngine::DenseMatrix cutset;
    SymEngine::Expression duty;
    std::vector<SymEngine::RCP<const SymEngine::Symbol>> symbols;
    SymEngine::DenseMatrix a_vector, b_vector, r_vector, ar_vector;

    SCCPhase(const SymEngine::DenseMatrix &on, const SymEngine::DenseMatrix &caps,
             const SymEngine::DenseMatrix &off, const SymEngine::DenseMatrix &loads,
             unsigned caps_count, const std::vector<unsigned> &indexes,
             const SymEngine::DenseMatrix &supply);
    SymEngine::DenseMatrix get_on_no_sw() const;
    void set_a_vector(const SymEngine::DenseMatrix &value) { a_vector = value; }
    const SymEngine::DenseMatrix &get_a_vector() const { return a_vector; }
};

#endif
