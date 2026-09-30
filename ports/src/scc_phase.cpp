#include "scc_phase.h"
#include <stdexcept>

using SymEngine::DenseMatrix;

static DenseMatrix concat(const DenseMatrix &a, const DenseMatrix &b) {
    if (a.nrows() != b.nrows()) throw std::invalid_argument("SCC_Phase: row dimension mismatch");
    DenseMatrix out(a.nrows(), a.ncols()+b.ncols());
    for (unsigned r=0;r<a.nrows();++r) { for (unsigned c=0;c<a.ncols();++c) out.set(r,c,a.get(r,c)); for (unsigned c=0;c<b.ncols();++c) out.set(r,a.ncols()+c,b.get(r,c)); }
    return out;
}

SCCPhase::SCCPhase(const DenseMatrix &on, const DenseMatrix &caps, const DenseMatrix &off,
                   const DenseMatrix &loads, unsigned caps_count,
                   const std::vector<unsigned> &indexes, const DenseMatrix &supply)
    : n_caps(caps_count), n_loads(loads.ncols()), n_on_sw(on.ncols()), n_off_sw(off.ncols()),
      inc_on_sw(on), sw_idxs(indexes) {
    if (supply.ncols()!=1) throw std::invalid_argument("SCC_Phase: supply must be one branch");
    inc_on_conv_sw = concat(concat(concat(supply,caps),on),loads);
    inc_on_conv = phase_conv(concat(concat(concat(supply,caps),loads),off),on);
    if (inc_on_conv.ncols() < n_off_sw) throw std::invalid_argument("SCC_Phase: invalid contracted dimensions");
}

DenseMatrix SCCPhase::get_on_no_sw() const {
    DenseMatrix out(inc_on_conv.nrows(), inc_on_conv.ncols()-n_off_sw);
    for (unsigned r=0;r<out.nrows();++r) for (unsigned c=0;c<out.ncols();++c) out.set(r,c,inc_on_conv.get(r,c));
    return out;
}
