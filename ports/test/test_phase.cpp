#include "dickson_hybrid_topology.h"
#include "phase_conv.h"
#include <cassert>
#include <symengine/integer.h>
int main() {
    Topology t = dickson_hybrid_topology(2, .5, true, false);
    assert(t.phase.size() == 2 && t.duty.ncols() == 2);
    assert(t.phase[0].sw_idxs == std::vector<unsigned>({1,3}));
    SymEngine::DenseMatrix e(2,2), s(2,1);
    e.set(0,0,SymEngine::integer(1)); e.set(1,0,SymEngine::integer(-1));
    e.set(0,1,SymEngine::integer(0)); e.set(1,1,SymEngine::integer(1));
    s.set(0,0,SymEngine::integer(1)); s.set(1,0,SymEngine::integer(-1));
    auto out=phase_conv(e,s); assert(out.nrows()==1 && out.ncols()==2);
    return 0;
}
