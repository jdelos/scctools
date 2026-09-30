#include "utilities.h"
#include <stdexcept>

SymEngine::DenseMatrix append_mA(const SymEngine::DenseMatrix &a, const SymEngine::DenseMatrix &b) {
    if (a.nrows() != b.nrows()) throw std::invalid_argument("append_mA: row mismatch");
    SymEngine::DenseMatrix out(a.nrows(), a.ncols()+b.ncols());
    for (unsigned r=0;r<a.nrows();++r) { for(unsigned c=0;c<a.ncols();++c) out.set(r,c,a.get(r,c)); for(unsigned c=0;c<b.ncols();++c) out.set(r,a.ncols()+c,b.get(r,c)); }
    return out;
}
