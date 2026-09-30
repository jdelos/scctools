#include "dickson_hybrid_topology.h"
#include <cassert>
#include <initializer_list>
#include <iostream>
#include <stdexcept>
#include <string>
#include <symengine/integer.h>
#include <symengine/parser.h>
#include <symengine/subs.h>
#include <symengine/symbol.h>

using SymEngine::DenseMatrix;
using SymEngine::Expression;
using SymEngine::RCP;
using SymEngine::Basic;

static DenseMatrix M(unsigned rows, unsigned cols, std::initializer_list<const char *> values) {
    assert(values.size() == rows * cols);
    DenseMatrix out(rows, cols); unsigned i = 0;
    for (const char *s : values) { out.set(i / cols, i % cols, SymEngine::parse(s)); ++i; }
    return out;
}
static void same(const DenseMatrix &a, const DenseMatrix &b) {
    assert(a.nrows() == b.nrows() && a.ncols() == b.ncols());
    for (unsigned r = 0; r < a.nrows(); ++r)
        for (unsigned c = 0; c < a.ncols(); ++c)
            if (!SymEngine::eq(*a.get(r,c), *b.get(r,c))) { std::cerr << r << "," << c << ": " << *a.get(r,c) << " != " << *b.get(r,c) << "\\n"; assert(false); }
}
static void expect_loads(const SCCPhase &p, const char *weight) {
    for (unsigned c = 0; c < p.n_loads; ++c)
        for (unsigned r = 0; r < p.inc_on_conv_sw.nrows(); ++r) {
            SymEngine::RCP<const SymEngine::Basic> expected = r == c + 1 ? SymEngine::parse(weight) : SymEngine::RCP<const SymEngine::Basic>(SymEngine::integer(0));
            assert(SymEngine::eq(*p.inc_on_conv_sw.get(r, p.inc_on_conv_sw.ncols() - p.n_loads + c), *expected));
        }
}
static void expect(const DenseMatrix &got, unsigned r, unsigned c,
                   std::initializer_list<const char *> v) { same(got, M(r,c,v)); }
static void expect_phase(const SCCPhase &p, unsigned n_caps, unsigned loads,
                         std::initializer_list<unsigned> indexes,
                         unsigned on, unsigned off, const DenseMatrix &sw,
                         const DenseMatrix &conv_sw, const DenseMatrix &conv,
                         const DenseMatrix &graph, const DenseMatrix &cutset) {
    assert(p.n_caps == n_caps && p.n_loads == loads && p.n_on_sw == on && p.n_off_sw == off);
    assert(p.sw_idxs == std::vector<unsigned>(indexes));
    same(p.inc_on_sw, sw); same(p.inc_on_conv_sw, conv_sw); same(p.inc_on_conv, conv);
    same(p.graph, graph); same(p.cutset, cutset);
}
static void rejected(int n, const Expression &d, std::vector<RCP<const SymEngine::Symbol>> s,
                     bool dc, bool half) {
    bool bad = false; try { (void)dickson_hybrid_topology(n,d,s,dc,half); }
    catch (const std::invalid_argument &) { bad = true; } assert(bad);
}

int main() {
    auto d = SymEngine::symbol("D"); Expression D(d);
    Topology t2 = dickson_hybrid_topology(2,D,{d},true,false);
    Topology t3 = dickson_hybrid_topology(3,D,{d},true,false);
    assert(t2.ordered_symbols == std::vector<RCP<const SymEngine::Symbol>>{d});
    assert(t3.ordered_symbols == t2.ordered_symbols);
    expect(t2.duty,1,2,{"D","1-D"}); expect(t3.duty,1,2,{"D","1-D"});

    DenseMatrix s2p0=M(4,2,{"1","0","-1","0","0","1","0","-1"});
    DenseMatrix s2p1=M(4,2,{"0","0","1","0","-1","0","0","1"});
    expect_phase(t2.phase[0],2,3,{1,3},2,2,s2p0,
      M(4,8,{"1","0","0","1","0","0","0","0","0","1","0","-1","0","D","0","0","0","0","1","0","1","0","D","0","0","-1","0","0","-1","0","0","D"}),
      M(2,8,{"1","1","0","D","0","0","1","0","0","-1","1","0","D","D","-1","1"}),
      M(2,6,{"1","1","0","D","0","0","0","-1","1","0","D","D"}),
      M(2,6,{"1","0","1","D","D","D","0","1","-1","0","-D","-D"}));
    expect_phase(t2.phase[1],2,3,{2,4},2,2,s2p1,
      M(4,8,{"1","0","0","0","0","0","0","0","0","1","0","1","0","1-D","0","0","0","0","1","-1","0","0","1-D","0","0","-1","0","0","1","0","0","1-D"}),
      M(2,8,{"1","0","0","0","0","0","1","0","0","1","1","1-D","1-D","0","-1","1"}),
      M(2,6,{"1","0","0","0","0","0","0","1","1","1-D","1-D","0"}),
      M(2,6,{"1","0","0","0","0","0","0","1","1","1-D","1-D","0"}));

    assert(t3.phase[0].sw_idxs == std::vector<unsigned>({1,3,5,7}));
    assert(t3.phase[1].sw_idxs == std::vector<unsigned>({2,4,6}));
    assert(t3.phase[0].graph.ncols() == 9 && t3.phase[1].graph.ncols() == 9);
    for (const Topology *tp : {&t2,&t3}) { assert(tp->phase.size()==2); for (const auto &p:tp->phase) { assert(p.tree.size()==p.graph.nrows()); assert(p.symbols.size()==1 && p.symbols[0]==d); } }
    expect_loads(t2.phase[0], "D"); expect_loads(t2.phase[1], "1-D");
    expect_loads(t3.phase[0], "D"); expect_loads(t3.phase[1], "1-D");

    Topology n = dickson_hybrid_topology(3,.25,true,false);
    assert(n.ordered_symbols.empty());
    assert(SymEngine::eq(*n.duty.get(0, 0), *SymEngine::real_double(.25)));
    assert(SymEngine::eq(*n.duty.get(0, 1), *SymEngine::real_double(.75)));
    for(unsigned p=0;p<2;++p) {
        assert(n.phase[p].graph.nrows() == t3.phase[p].graph.nrows());
        assert(n.phase[p].graph.ncols() == t3.phase[p].graph.ncols());
        assert(n.phase[p].cutset.nrows() == t3.phase[p].cutset.nrows());
        assert(n.phase[p].cutset.ncols() == t3.phase[p].cutset.ncols());
    }
    expect_loads(n.phase[0], "0.25"); expect_loads(n.phase[1], "0.75");
    rejected(2,D,{d,SymEngine::symbol("E")},true,false);
    rejected(1,D,{d},true,false); rejected(2,D,{},true,false); rejected(2,D,{d},false,false); rejected(2,D,{d},true,true);
    return 0;
}
