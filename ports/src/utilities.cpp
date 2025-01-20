#include "utilities.h"

// TODO: Implement function logic here

SymEngine::DenseMatrix append_mA(const SymEngine::DenseMatrix &A1, const SymEngine::DenseMatrix &A2) {
    SymEngine::DenseMatrix result(A1.nrows() + A2.nrows() - 1, A1.ncols() + A2.ncols());
    result.submatrix(0, 0, A1);
    result.submatrix(A1.nrows() - 1, A1.ncols(), A2);
    return result;
}