#ifndef SOLVE_CHARGE_VECTORS_H
#define SOLVE_CHARGE_VECTORS_H

#include <symengine/matrix.h>
#include <symengine/expression.h>
#include <symengine/symbol.h>
#include <vector>

struct ChargeSolution {
    SymEngine::DenseMatrix Qx;
    SymEngine::DenseMatrix Qo;
    std::vector<SymEngine::DenseMatrix> a;
    SymEngine::DenseMatrix m;
};

ChargeSolution solve_charge_vectors(
    const std::vector<SymEngine::DenseMatrix> &cutsets,
    unsigned n_caps,
    const SymEngine::DenseMatrix &duty,
    const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &symbols);

#endif
