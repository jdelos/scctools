#include "dickson_hybrid_topology.h"
#include "graph_primitives.h"
#include "utilities.h"
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <stdexcept>
using SymEngine::DenseMatrix;

Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point) {
    ArchDef arch = dickson_arch(n_caps);

    if (half_point) {
        SymEngine::DenseMatrix A_cap_hp(4, 1), A_sw1_hp(4, 2), A_sw2_hp(4, 2);
        append_mA(A_cap_hp, arch.Acaps);
        append_mA(A_sw1_hp, arch.Asw);
        append_mA(A_sw2_hp, arch.Asw);
    }

    if (dc_out == false || half_point) throw std::invalid_argument("dickson_hybrid_topology: unsupported option");
    Topology top;
    top.N_caps = n_caps; top.N_sw = arch.Asw.ncols(); top.vo_swing = 1.0 / n_caps;
    top.duty = SymEngine::DenseMatrix(1, 2);
    auto d = SymEngine::symbol("D");
    top.duty.set(0, 0, d); top.duty.set(0, 1, SymEngine::sub(SymEngine::integer(1), d));
    DenseMatrix loads(arch.Acaps.nrows(), arch.Acaps.nrows()-1), supply(arch.Acaps.nrows(),1);
    for (unsigned r=0;r<loads.nrows();++r) for (unsigned c=0;c<loads.ncols();++c) loads.set(r,c,SymEngine::integer(r==c+1));
    for (unsigned r=0;r<supply.nrows();++r) supply.set(r,0,SymEngine::integer(r==0));
    for (unsigned p=0;p<2;++p) {
        DenseMatrix on(arch.Asw.nrows(),0), off(arch.Asw.nrows(),0); std::vector<unsigned> indexes;
        for (unsigned c=0;c<arch.Asw.ncols();++c) {
            DenseMatrix one(arch.Asw.nrows(),1); for (unsigned r=0;r<one.nrows();++r) one.set(r,0,arch.Asw.get(r,c));
            bool active = SymEngine::eq(*arch.Asw_act.get(p,c), *SymEngine::integer(1));
            if (active) { DenseMatrix x(on.nrows(),on.ncols()+1); for(unsigned r=0;r<x.nrows();++r){for(unsigned j=0;j<on.ncols();++j)x.set(r,j,on.get(r,j));x.set(r,on.ncols(),one.get(r,0));} on=x; indexes.push_back(c+1); }
            else { DenseMatrix x(off.nrows(),off.ncols()+1); for(unsigned r=0;r<x.nrows();++r){for(unsigned j=0;j<off.ncols();++j)x.set(r,j,off.get(r,j));x.set(r,off.ncols(),one.get(r,0));} off=x; }
        }
        top.phase.emplace_back(on,arch.Acaps,off,loads,n_caps,indexes,supply);
        auto graph=top.phase.back().get_on_no_sw();
        top.phase.back().tree=build_tree(graph,0);
        if (top.phase.back().tree.size()==1 && top.phase.back().tree[0]==static_cast<unsigned>(-1)) throw std::runtime_error("dickson_hybrid_topology: singular graph");
        top.phase.back().cutset=fun_cutset(graph,top.phase.back().tree);
    }
    return top;
}
