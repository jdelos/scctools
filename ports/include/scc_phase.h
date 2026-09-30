#ifndef SCC_PHASE_H
#define SCC_PHASE_H

#include "graph_primitives.h"
#include "phase_conv.h"
#include <symengine/symbol.h>
#include <vector>

struct SCCPhase {
    unsigned n_caps, n_loads, n_on_sw, n_off_sw;
    SymEngine::DenseMatrix inc_on_sw, inc_on_conv, inc_on_conv_sw;
    std::vector<unsigned> sw_idxs;
    std::vector<unsigned> tree;
    SymEngine::DenseMatrix cutset;

    SCCPhase(const SymEngine::DenseMatrix &on, const SymEngine::DenseMatrix &caps,
             const SymEngine::DenseMatrix &off, const SymEngine::DenseMatrix &loads,
             unsigned caps_count, const std::vector<unsigned> &indexes,
             const SymEngine::DenseMatrix &supply);
    SymEngine::DenseMatrix get_on_no_sw() const;
};

#endif
