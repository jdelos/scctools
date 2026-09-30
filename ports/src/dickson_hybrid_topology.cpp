#include "dickson_hybrid_topology.h"
#include "graph_primitives.h"
#include "utilities.h"
#include "solve_charge_vectors.h"
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/visitor.h>
#include <stdexcept>
using SymEngine::DenseMatrix;

Topology dickson_hybrid_topology(int n_caps, const SymEngine::Expression &duty,
                                 const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &symbols,
                                 bool dc_out, bool half_point) {
    if (dc_out == false || half_point) throw std::invalid_argument("dickson_hybrid_topology: unsupported option");
    auto free = SymEngine::free_symbols(*duty.get_basic());
    if (free.size() != symbols.size()) throw std::invalid_argument("dickson_hybrid_topology: symbols do not match duty");
    for (const auto &symbol : symbols)
        if (free.find(symbol) == free.end()) throw std::invalid_argument("dickson_hybrid_topology: symbols do not match duty");
    ArchDef arch = dickson_arch(n_caps);
    Topology top;
    top.ordered_symbols = symbols;
    top.N_caps = n_caps; top.N_sw = arch.Asw.ncols(); top.vo_swing = 1.0 / n_caps;
    top.duty = SymEngine::DenseMatrix(1, 2);
    top.duty.set(0, 0, duty.get_basic());
    top.duty.set(0, 1, SymEngine::sub(SymEngine::integer(1), duty.get_basic()));
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
        DenseMatrix phase_loads = loads;
        auto phase_duty = p == 0 ? duty.get_basic() : SymEngine::sub(SymEngine::integer(1), duty.get_basic());
        for (unsigned c = 0; c < phase_loads.ncols(); ++c)
            phase_loads.set(c + 1, c, phase_duty);
        top.phase.emplace_back(on,arch.Acaps,off,phase_loads,n_caps,indexes,supply);
        top.phase.back().duty = p == 0 ? duty : SymEngine::Expression(SymEngine::sub(SymEngine::integer(1), duty.get_basic()));
        top.phase.back().symbols = symbols;
        top.phase.back().graph = top.phase.back().get_on_no_sw();
        auto graph=top.phase.back().graph;
        top.phase.back().tree=build_tree(graph,0);
        if (top.phase.back().tree.size()==1 && top.phase.back().tree[0]==static_cast<unsigned>(-1)) throw std::runtime_error("dickson_hybrid_topology: singular graph");
        top.phase.back().cutset=fun_cutset(graph,top.phase.back().tree);
    }
    std::vector<DenseMatrix> cutsets;
    for (const auto &phase : top.phase) cutsets.push_back(phase.cutset);
    ChargeSolution charge = solve_charge_vectors(cutsets, n_caps, top.duty, symbols);
    top.m_ratios = charge.m;
    for (unsigned p = 0; p < top.phase.size(); ++p) top.phase[p].set_a_vector(charge.a[p]);
    return top;
}

Topology dickson_hybrid_topology(int n_caps, double duty, bool dc_out, bool half_point) {
    return dickson_hybrid_topology(n_caps, SymEngine::Expression(duty),
                                   {}, dc_out, half_point);
}
