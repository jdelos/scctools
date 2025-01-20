#include "dickson_hybrid_topology.h"

Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point) {
    ArchDef arch = dickson_arch(n_caps);

    if (half_point) {
        SymEngine::DenseMatrix A_cap_hp(4, 1), A_sw1_hp(4, 2), A_sw2_hp(4, 2);
        append_mA(A_cap_hp, arch.Acaps);
        append_mA(A_sw1_hp, arch.Asw);
        append_mA(A_sw2_hp, arch.Asw);
    }

    Topology top;
    top.N_caps = n_caps;
    top.duty = duty;
    top.vo_swing = 1.0 / n_caps;

    return top;
}
