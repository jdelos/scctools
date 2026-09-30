#include "graph_primitives.h"
#include <symengine/integer.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <algorithm>
#include <functional>
#include <stdexcept>

using SymEngine::Basic;
using SymEngine::DenseMatrix;
using SymEngine::RCP;

static bool nonzero(const RCP<const Basic> &x) {
    return !SymEngine::eq(*x, *SymEngine::integer(0));
}

static DenseMatrix with_ground(const DenseMatrix &a) {
    DenseMatrix out(a.nrows() + 1, a.ncols());
    for (unsigned r = 0; r < a.nrows(); ++r)
        for (unsigned c = 0; c < a.ncols(); ++c) out.set(r, c, a.get(r, c));
    for (unsigned c = 0; c < a.ncols(); ++c) {
        RCP<const Basic> sum = SymEngine::integer(0);
        for (unsigned r = 0; r < a.nrows(); ++r) sum = SymEngine::add(sum, a.get(r, c));
        out.set(a.nrows(), c, SymEngine::neg(sum));
    }
    return out;
}

static bool parallel(const DenseMatrix &a, unsigned x, unsigned y) {
    for (unsigned r = 0; r < a.nrows(); ++r)
        if (nonzero(a.get(r, x)) != nonzero(a.get(r, y))) return false;
    // Incidence columns are compared by absolute value in MATLAB. For the
    // numeric/symbolic values used here, equality of squares avoids sign.
    for (unsigned r = 0; r < a.nrows(); ++r) {
        auto xx = a.get(r, x), yy = a.get(r, y);
        if (!SymEngine::eq(*SymEngine::mul(xx, xx), *SymEngine::mul(yy, yy))) return false;
    }
    return true;
}

std::vector<unsigned> build_tree(const DenseMatrix &incidence, int initial_edge,
                                 const std::vector<unsigned> &excluded) {
    if (!incidence.nrows() || !incidence.ncols())
        throw std::invalid_argument("build_tree: empty incidence matrix");
    std::vector<bool> omit(incidence.ncols(), false);
    for (unsigned e : excluded) {
        if (e >= incidence.ncols()) throw std::invalid_argument("build_tree: excluded index");
        omit[e] = true;
    }
    DenseMatrix a = with_ground(incidence);
    std::vector<unsigned> columns;
    for (unsigned e = 0; e < incidence.ncols(); ++e) if (!omit[e]) columns.push_back(e);
    unsigned start;
    if (initial_edge < 0) {
        auto it = std::find_if(columns.begin(), columns.end(), [&](unsigned e) {
            for (unsigned r = 0; r < a.nrows(); ++r) if (nonzero(a.get(r, e))) return true;
            return false;
        });
        if (it == columns.end()) throw std::runtime_error("build_tree: graph has no spanning tree");
        start = *it;
    } else {
        start = static_cast<unsigned>(initial_edge);
        if (start >= incidence.ncols() || omit[start])
            throw std::invalid_argument("build_tree: invalid initial edge");
    }
    bool valid = false;
    for (unsigned r = 0; r < a.nrows(); ++r) valid |= nonzero(a.get(r, start));
    if (!valid) throw std::invalid_argument("build_tree: initial edge is not valid");

    const unsigned need = a.nrows() - 1;
    std::vector<unsigned> tree{start};
    std::vector<bool> listed(a.nrows(), false);
    for (unsigned r = 0; r < a.nrows(); ++r) listed[r] = nonzero(a.get(r, start));

    std::function<bool(unsigned, std::vector<unsigned>, std::vector<bool>, DenseMatrix)> visit;
    visit = [&](unsigned incoming, std::vector<unsigned> path, std::vector<bool> nodes,
                DenseMatrix state) -> bool {
        if (path.size() == need) { tree = std::move(path); return true; }
        unsigned current = path.back();
        for (unsigned e : columns)
            if (e != current && parallel(state, current, e))
                for (unsigned r = 0; r < state.nrows(); ++r) state.set(r, e, SymEngine::integer(0));
        for (unsigned r = 0; r < state.nrows(); ++r) state.set(r, current, SymEngine::integer(0));
        std::vector<unsigned> ends;
        for (unsigned r = 0; r < state.nrows(); ++r)
            if (nonzero(a.get(r, current)) && r != incoming) ends.push_back(r);
        for (unsigned node : ends) {
            for (unsigned e : columns) {
                if (nonzero(state.get(node, e))) {
                    unsigned hits = 0;
                    for (unsigned r = 0; r < state.nrows(); ++r) hits += nonzero(state.get(r, e));
                    if (hits != 2 || std::find(path.begin(), path.end(), e) != path.end()) continue;
                    unsigned other_hits = 0;
                    for (unsigned r = 0; r < state.nrows(); ++r)
                        if (nonzero(a.get(r, e)) && nodes[r]) ++other_hits;
                    if (other_hits == 2) continue;
                    auto next = path; next.push_back(e);
                    auto next_nodes = nodes;
                    for (unsigned r = 0; r < state.nrows(); ++r) if (nonzero(a.get(r, e))) next_nodes[r] = true;
                    if (visit(node, std::move(next), std::move(next_nodes), state)) return true;
                }
            }
        }
        return false;
    };
    if (!visit(a.nrows(), tree, listed, a)) return {-1u};
    std::sort(tree.begin(), tree.end());
    return tree;
}

bool full_tree(const DenseMatrix &a) {
    if (!a.nrows() || !a.ncols()) return false;
    // MATLAB full_tree checks rows for matrices, but entries for column vectors.
    if (a.ncols() == 1) {
        for (unsigned r = 0; r < a.nrows(); ++r) if (!nonzero(a.get(r, 0))) return false;
    } else {
        for (unsigned r = 0; r < a.nrows(); ++r) {
            bool present = false;
            for (unsigned c = 0; c < a.ncols(); ++c) present |= nonzero(a.get(r, c));
            if (!present) return false;
        }
    }
    return true;
}

DenseMatrix fun_cutset(const DenseMatrix &a, const std::vector<unsigned> &requested) {
    if (!a.nrows() || !a.ncols()) throw std::invalid_argument("fun_cutset: empty matrix");
    std::vector<unsigned> tree = requested.empty() ? build_tree(a, 0) : requested;
    if (tree.size() != a.nrows()) throw std::invalid_argument("fun_cutset: tree dimension");
    std::vector<bool> used(a.ncols(), false);
    for (unsigned e : tree) {
        if (e >= a.ncols() || used[e]) throw std::invalid_argument("fun_cutset: invalid tree index");
        used[e] = true;
    }
    std::vector<unsigned> order = tree;
    for (unsigned e = 0; e < a.ncols(); ++e) if (!used[e]) order.push_back(e);
    DenseMatrix at(a.nrows(), tree.size()), al(a.nrows(), a.ncols() - tree.size());
    for (unsigned r = 0; r < a.nrows(); ++r) {
        for (unsigned c = 0; c < order.size(); ++c)
            (c < tree.size() ? at : al).set(r, c < tree.size() ? c : c - tree.size(), a.get(r, order[c]));
    }
    DenseMatrix x(tree.size(), al.ncols());
    at.LU_solve(al, x);
    DenseMatrix q_ordered(a.nrows(), a.ncols());
    for (unsigned r = 0; r < a.nrows(); ++r) {
        for (unsigned c = 0; c < tree.size(); ++c) q_ordered.set(r, c, SymEngine::integer(r == c));
        for (unsigned c = tree.size(); c < order.size(); ++c) q_ordered.set(r, c, x.get(r, c - tree.size()));
    }
    DenseMatrix q(a.nrows(), a.ncols());
    for (unsigned c = 0; c < order.size(); ++c) for (unsigned r = 0; r < a.nrows(); ++r) q.set(r, order[c], q_ordered.get(r, c));
    return q;
}
