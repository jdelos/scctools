#include "dickson_hybrid_topology.h"
#include <stdexcept>

Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point) {
    if (n_caps < 2) {
        throw std::invalid_argument("invalid architecture: n_caps must be at least 2");
    }
    if (duty <= 0.0 || duty >= 1.0) {
        throw std::invalid_argument("singular input: duty must be between zero and one");
    }
    if (half_point) {
        throw std::invalid_argument("half-point topology unsupported at two-phase boundary");
    }
    (void)dc_out;
    ArchDef arch = dickson_arch(n_caps);
    Topology top;
    top.A = arch.A;
    top.m = arch.m;
    top.ratio = arch.m;
    top.N_caps = n_caps;
    top.N_sw = arch.Asw.ncols();
    top.duty = duty;
    top.vo_swing = 1.0 / n_caps;
    return top;
}
