#ifndef NATIVE_GRAPH_H
#define NATIVE_GRAPH_H
#include <symengine/matrix.h>
#include <string>
#include <vector>

struct NativePhase {
    SymEngine::DenseMatrix graph;
    SymEngine::DenseMatrix cutset;
    SymEngine::DenseMatrix inc_on_conv;
    SymEngine::DenseMatrix inc_on_conv_sw;
    std::vector<unsigned> tree;
    std::vector<unsigned> switch_indices;
    unsigned n_on_sw = 0;
    unsigned n_off_sw = 0;
};
struct NativeGraph {
    std::vector<NativePhase> phases;
    std::vector<std::string> ordered_symbols;
    SymEngine::DenseMatrix duties;
};
// Native graph slice: exactly two phases, symbolic duty D and 1-D.
NativeGraph native_graph(int n_caps);
#endif
