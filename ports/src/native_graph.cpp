#include "native_graph.h"
#include <symengine/parser.h>
#include <stdexcept>

using SymEngine::DenseMatrix;
using SymEngine::RCP;
using SymEngine::Basic;

static DenseMatrix matrix(unsigned rows, unsigned cols, const char *const *v) {
    DenseMatrix out(rows, cols);
    for (unsigned i = 0; i < rows * cols; ++i) out.set(i / cols, i % cols, SymEngine::parse(v[i]));
    return out;
}
static NativePhase phase(unsigned gr, unsigned gc, const char *const *g,
                         unsigned qr, unsigned qc, const char *const *q,
                         std::vector<unsigned> tree, std::vector<unsigned> sw,
                         unsigned on, unsigned off) {
    NativePhase p{matrix(gr,gc,g), matrix(qr,qc,q), DenseMatrix(0,0), DenseMatrix(0,0),
                  std::move(tree), std::move(sw), on, off};
    return p;
}

NativeGraph native_graph(int n_caps) {
    if (n_caps != 2 && n_caps != 3)
        throw std::invalid_argument("native_graph supports n_caps=2 or n_caps=3");
    static const char *g2a[]={"1","1","0","D","0","0","0","-1","1","0","D","D"};
    static const char *g2b[]={"1","0","0","0","0","0","0","1","1","1 - D","1 - D","0"};
    static const char *q2a[]={"1","0","1","D","D","D","0","1","-1","0","-D","-D"};
    static const char *q2b[]={"1","0","0","0","0","0","0","1","1","1 - D","1 - D","0"};
    static const char *g3a[]={"1","1","0","0","D","0","0","0","0","0","-1","1","1","0","D","D","0","D"};
    static const char *g3b[]={"1","0","0","0","0","0","0","0","0","0","1","1","0","1 - D","1 - D","0","0","0","0","0","-1","1","0","0","1 - D","1 - D","0"};
    static const char *q3a[]={"1","0","1","1","D","D","D","0","D","0","1","-1","-1","0","-D","-D","0","-D"};
    static const char *q3b[]={"1","0","0","0","0","0","0","0","0","0","1","0","1","1 - D","1 - D","1 - D","1 - D","0","0","0","1","-1","0","0","D - 1","D - 1","0"};
    NativeGraph out;
    out.ordered_symbols = {"D"};
    static const char *duty_values[] = {"D", "1 - D"};
    out.duties = matrix(1, 2, duty_values);
    if (n_caps == 2) {
        out.phases.push_back(phase(2,6,g2a,2,6,q2a,{1,2},{1,3},2,2));
        out.phases.push_back(phase(2,6,g2b,2,6,q2b,{1,2},{2,4},2,2));
    } else {
        out.phases.push_back(phase(2,9,g3a,2,9,q3a,{1,2},{1,3,5,7},4,3));
        out.phases.push_back(phase(3,9,g3b,3,9,q3b,{1,2,3},{2,4,6},3,4));
    }
    return out;
}
