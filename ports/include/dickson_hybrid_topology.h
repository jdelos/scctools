
#ifndef DICKSON_HYBRID_TOPOLOGY_H
#define DICKSON_HYBRID_TOPOLOGY_H

#include "dickson_arch.h"
#include <symengine/expression.h>
#include <vector>
#include "scc_phase.h"

struct Topology {
    SymEngine::DenseMatrix ratio, vc, vr, is, m_ratios;
    SymEngine::DenseMatrix ZSSL, ZFSL, ZESR, ZSCC;
    std::vector<SymEngine::RCP<const SymEngine::Symbol>> capacitances, switch_resistances, capacitor_esr;
    SymEngine::RCP<const SymEngine::Symbol> frequency;
    double vo_swing;
    int N_caps, N_sw;
    SymEngine::DenseMatrix duty;
    std::vector<SCCPhase> phase;
    std::vector<SymEngine::RCP<const SymEngine::Symbol>> ordered_symbols;
};

Topology dickson_hybrid_topology(int n_caps, const SymEngine::Expression &duty,
                                 const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &symbols,
                                 bool dc_out, bool half_point);
Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point);

#endif

