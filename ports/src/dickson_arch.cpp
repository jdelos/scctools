#include "dickson_arch.h"
#include "dickson_matrix.h"
#include <stdexcept>
#include <vector>

using SymEngine::DenseMatrix;
using SymEngine::integer;

namespace {
DenseMatrix to_symbolic(const std::vector<std::vector<int>> &source) {
    DenseMatrix result(source.size(), source.empty() ? 0 : source.front().size());
    for (unsigned row = 0; row < source.size(); ++row)
        for (unsigned col = 0; col < source[row].size(); ++col)
            result.set(row, col, integer(source[row][col]));
    return result;
}
}

ArchDef dickson_arch(int n_caps) {
    if (n_caps < 2) {
        throw std::invalid_argument("invalid architecture: n_caps must be at least 2");
    }
    std::vector<std::vector<int>> caps, phase1, phase2;
    dickson_matrix(n_caps, false, caps, phase1, phase2);

    ArchDef result;
    result.Acaps = to_symbolic(caps);
    const unsigned phase1_cols = phase1.front().size();
    const unsigned phase2_cols = phase2.front().size();
    result.Asw = DenseMatrix(caps.size(), phase1_cols + phase2_cols);
    result.Asw_act = DenseMatrix(2, phase1_cols + phase2_cols);
    for (unsigned row = 0; row < caps.size(); ++row) {
        for (unsigned col = 0; col < phase1_cols; ++col)
            result.Asw.set(row, col, integer(phase1[row][col]));
        for (unsigned col = 0; col < phase2_cols; ++col)
            result.Asw.set(row, phase1_cols + col, integer(phase2[row][col]));
    }
    for (unsigned col = 0; col < phase1_cols; ++col)
        result.Asw_act.set(0, col, integer(1));
    for (unsigned col = 0; col < phase2_cols; ++col)
        result.Asw_act.set(1, phase1_cols + col, integer(1));
    result.A = result.Acaps;
    result.m = DenseMatrix(1, n_caps - 1);
    for (unsigned col = 0; col < result.m.ncols(); ++col)
        result.m.set(0, col, integer(static_cast<int>(col + 2)));
    return result;
}
