#ifndef GRAPH_PRIMITIVES_H
#define GRAPH_PRIMITIVES_H

#include <symengine/matrix.h>
#include <vector>

// Native graph indexes are zero-based. Excluded indexes refer to original columns.
std::vector<unsigned> build_tree(const SymEngine::DenseMatrix &incidence,
                                 int initial_edge = -1,
                                 const std::vector<unsigned> &excluded = {});
bool full_tree(const SymEngine::DenseMatrix &tree_incidence);
SymEngine::DenseMatrix fun_cutset(const SymEngine::DenseMatrix &incidence,
                                  const std::vector<unsigned> &tree = {});

#endif
