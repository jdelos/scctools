#include "dickson_hybrid_topology.h"
#include "phase_conv.h"
#include <cassert>
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/rational.h>
#include <stdexcept>

static void assert_matrix(const SymEngine::DenseMatrix &a, const SymEngine::DenseMatrix &b) {
    assert(a.nrows() == b.nrows() && a.ncols() == b.ncols());
    for (unsigned r = 0; r < a.nrows(); ++r)
        for (unsigned c = 0; c < a.ncols(); ++c)
            assert(SymEngine::eq(*a.get(r, c), *b.get(r, c)));
}

int main() {
    auto d = SymEngine::symbol("D");
    Topology t = dickson_hybrid_topology(2, SymEngine::Expression(d), {d}, true, false);
    assert(t.phase.size() == 2 && t.duty.nrows() == 1 && t.duty.ncols() == 2);
    assert(SymEngine::eq(*t.duty.get(0, 0), *d));
    assert(SymEngine::eq(*t.duty.get(0, 1), *SymEngine::sub(SymEngine::integer(1), d)));
    assert(t.phase[0].sw_idxs == std::vector<unsigned>({1,3}));
    assert(t.phase[1].sw_idxs == std::vector<unsigned>({2,4}));
    assert(t.phase[0].n_caps == 2 && t.phase[0].n_loads == 3);
    assert(t.phase[0].n_on_sw == 2 && t.phase[0].n_off_sw == 2);
    assert(t.phase[0].symbols.size() == 1 && SymEngine::eq(*t.phase[0].symbols[0], *d));
    assert(SymEngine::eq(*t.phase[0].duty.get_basic(), *d));
    assert(SymEngine::eq(*t.phase[1].duty.get_basic(), *SymEngine::sub(SymEngine::integer(1), d)));
    assert_matrix(t.phase[0].graph, t.phase[0].get_on_no_sw());
    assert_matrix(t.phase[1].graph, t.phase[1].get_on_no_sw());
    assert(t.phase[0].tree.size() == t.phase[1].tree.size());
    assert(t.phase[0].tree.size() == 2);
    assert(t.phase[0].tree == std::vector<unsigned>({0,1}));
    assert(t.phase[1].tree == std::vector<unsigned>({0,1}));
    assert(t.phase[0].cutset.nrows() == t.phase[0].tree.size());
    assert(t.phase[0].cutset.ncols() == t.phase[0].graph.ncols());
    SymEngine::DenseMatrix av(3, 1); av.set(0, 0, SymEngine::integer(1));
    t.phase[0].set_a_vector(av); assert_matrix(t.phase[0].get_a_vector(), av);
    bool rejected = false;
    try { (void)dickson_hybrid_topology(2, SymEngine::Expression(d), {d}, false, false); }
    catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)dickson_hybrid_topology(2, SymEngine::Expression(d), {d}, true, true); }
    catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    Topology numeric = dickson_hybrid_topology(2, .25, true, false);
    assert(SymEngine::eq(*numeric.duty.get(0, 0), *SymEngine::real_double(.25)));
    SymEngine::DenseMatrix e(2,2), s(2,1);
    e.set(0,0,SymEngine::integer(1)); e.set(1,0,SymEngine::integer(-1));
    e.set(0,1,SymEngine::integer(0)); e.set(1,1,SymEngine::integer(1));
    s.set(0,0,SymEngine::integer(1)); s.set(1,0,SymEngine::integer(-1));
    auto out=phase_conv(e,s); assert(out.nrows()==1 && out.ncols()==2);
    return 0;
}
