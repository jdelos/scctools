#include "graph_primitives.h"
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/add.h>
#include <cassert>
#include <stdexcept>

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
    DenseMatrix signed_cutset(2, 3);
    signed_cutset.set(0, 0, integer(1)); signed_cutset.set(1, 0, integer(-1));
    signed_cutset.set(0, 1, integer(-1)); signed_cutset.set(1, 1, integer(0));
    signed_cutset.set(0, 2, integer(0)); signed_cutset.set(1, 2, integer(1));
    const DenseMatrix signed_q = fun_cutset(signed_cutset, {1, 2});
    assert(eq(*signed_q.get(0, 0), *integer(-1)));
    bool cutset_rejected = false;
    try { (void)fun_cutset(a, {0, 0}); } catch (const std::invalid_argument &) { cutset_rejected = true; }
    assert(cutset_rejected);

    DenseMatrix symbolic(2, 3);
    auto x = symbol("x");
    symbolic.set(0, 0, integer(1)); symbolic.set(1, 0, integer(-1));
    symbolic.set(0, 1, add(x, integer(-1))); symbolic.set(1, 1, integer(0));
    symbolic.set(0, 2, integer(1)); symbolic.set(1, 2, integer(0));
    // x - 1 becomes zero in topology view; default selection skips edge 1.
    assert((build_tree(symbolic) == std::vector<unsigned>{0, 2}));
    DenseMatrix parallel_symbolic(2, 3);
    parallel_symbolic.set(0, 0, integer(1)); parallel_symbolic.set(1, 0, integer(-1));
    parallel_symbolic.set(0, 1, x); parallel_symbolic.set(1, 1, integer(-1));
    parallel_symbolic.set(0, 2, integer(0)); parallel_symbolic.set(1, 2, integer(1));
    assert((build_tree(parallel_symbolic) == std::vector<unsigned>{0, 2}));
    assert(full_tree(symbolic));
    assert((build_tree(a, 0, {2}) == std::vector<unsigned>{0, 1}));
    bool rejected = false;
    try { (void)build_tree(a, 0, {9}); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)build_tree(a, 9); } catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);

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
