#include "graph_primitives.h"
#include <symengine/integer.h>
#include <cassert>

using namespace SymEngine;

int main() {
    DenseMatrix a(2, 3);
    a.set(0, 0, integer(1)); a.set(1, 0, integer(-1));
    a.set(0, 1, integer(0)); a.set(1, 1, integer(1));
    a.set(0, 2, integer(1)); a.set(1, 2, integer(0));
    const std::vector<unsigned> tree = build_tree(a, 0);
    assert((tree == std::vector<unsigned>{0, 2}));
    assert((build_tree(a) == std::vector<unsigned>{0, 2}));
    assert((build_tree(a, 0, {1}) == std::vector<unsigned>{0, 2}));
    assert(full_tree(a));
    const DenseMatrix q = fun_cutset(a, tree);
    assert(eq(*q.get(0, 0), *integer(1)));
    assert(eq(*q.get(1, 1), *integer(1)));
    assert(eq(*q.get(0, 2), *integer(0)));
    assert(eq(*q.get(1, 2), *integer(1)));

    DenseMatrix symbolic(2, 2);
    symbolic.set(0, 0, integer(1)); symbolic.set(1, 0, integer(-1));
    symbolic.set(0, 1, integer(1)); symbolic.set(1, 1, integer(0));
    assert(full_tree(symbolic));

    DenseMatrix disconnected(2, 2);
    for (unsigned r = 0; r < 2; ++r)
        for (unsigned c = 0; c < 2; ++c) disconnected.set(r, c, integer(0));
    disconnected.set(0, 0, integer(1)); disconnected.set(1, 0, integer(-1));
    assert(build_tree(disconnected) == std::vector<unsigned>{-1u});

    DenseMatrix backtrack(4, 5);
    int values[] = {1, 0, -1, 0, 0, -1, 1, 0, 0, 0, 0, -1, 1, 1, 0, 0, 0, 0, -1, 1};
    for (unsigned r = 0; r < 4; ++r)
        for (unsigned c = 0; c < 5; ++c) backtrack.set(r, c, integer(values[r * 5 + c]));
    assert((build_tree(backtrack) == std::vector<unsigned>{0, 2, 3, 4}));
    return 0;
}
