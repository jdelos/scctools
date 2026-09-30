
#ifndef DICKSON_HYBRID_TOPOLOGY_H
#define DICKSON_HYBRID_TOPOLOGY_H

#include "dickson_arch.h"
#include <symengine/expression.h>

struct Topology {
    SymEngine::DenseMatrix ratio, vc, vr, is;
    double vo_swing, duty;
    int N_caps, N_sw;
};

Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point);

#endif

