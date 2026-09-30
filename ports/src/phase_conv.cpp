#include "phase_conv.h"
#include <symengine/integer.h>
#include <symengine/add.h>
#include <stdexcept>
#include <vector>

using SymEngine::DenseMatrix;
using SymEngine::integer;

DenseMatrix phase_conv(const DenseMatrix &Aelem, const DenseMatrix &Asw) {
    if (Asw.nrows() != Aelem.nrows()) throw std::invalid_argument("phase_conv: row dimension mismatch");
    const unsigned m = Asw.ncols();
    DenseMatrix work(Asw.nrows(), m + Aelem.ncols());
    for (unsigned r = 0; r < work.nrows(); ++r) {
        for (unsigned c = 0; c < m; ++c) work.set(r, c, Asw.get(r, c));
        for (unsigned c = 0; c < Aelem.ncols(); ++c) work.set(r, m + c, Aelem.get(r, c));
    }
    unsigned rows = work.nrows();
    for (unsigned c = 0; c < m; ++c) {
        int pos = -1, neg = -1;
        for (unsigned r = 0; r < rows; ++r) {
            if (SymEngine::eq(*work.get(r, c), *integer(1))) pos = static_cast<int>(r);
            if (SymEngine::eq(*work.get(r, c), *integer(-1))) neg = static_cast<int>(r);
        }
        if (pos < 0 && neg < 0) throw std::invalid_argument("phase_conv: switch has no endpoint");
        unsigned remove = pos < 0 ? static_cast<unsigned>(neg) : static_cast<unsigned>(pos);
        if (pos >= 0 && neg >= 0) {
            for (unsigned j = 0; j < work.ncols(); ++j)
                work.set(static_cast<unsigned>(pos), j, SymEngine::add(work.get(pos,j), work.get(neg,j)));
            remove = static_cast<unsigned>(neg);
        }
        DenseMatrix next(rows - 1, work.ncols());
        for (unsigned r = 0, nr = 0; r < rows; ++r) if (r != remove) {
            for (unsigned j = 0; j < work.ncols(); ++j) next.set(nr,j,work.get(r,j));
            ++nr;
        }
        work = next; --rows;
    }
    DenseMatrix out(rows, Aelem.ncols());
    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < Aelem.ncols(); ++c) out.set(r,c,work.get(r,m+c));
    return out;
}
