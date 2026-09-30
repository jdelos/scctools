#include "native_graph.h"
#include <cassert>
#include <stdexcept>
#include <symengine/basic.h>
#include <symengine/integer.h>

static void check(int n, unsigned rows, unsigned cols) {
    NativeGraph g = native_graph(n);
    assert(g.phases.size() == 2 && g.ordered_symbols == std::vector<std::string>{"D"});
    assert(g.duties.nrows() == 1 && g.duties.ncols() == 2);
    assert(g.phases[0].graph.nrows() == rows && g.phases[0].graph.ncols() == cols);
    assert(g.phases[0].cutset.nrows() == rows && g.phases[0].cutset.ncols() == cols);
    assert(g.phases[0].graph.get(0, 3)->__str__().size() > 0);
}
int main() {
    check(2, 2, 6); check(3, 2, 9);
    bool rejected = false;
    try { native_graph(4); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
}
