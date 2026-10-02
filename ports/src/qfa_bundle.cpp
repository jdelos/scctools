#include "qfa_bundle.h"
#include <stdexcept>
#include <symengine/integer.h>
#include <symengine/pow.h>
#include <symengine/simplify.h>
#include <symengine/symbol.h>

using namespace SymEngine;
namespace {
DenseMatrix zeros(unsigned r, unsigned c) {
  DenseMatrix m(r, c);
  for (unsigned i = 0; i < r; ++i)
    for (unsigned j = 0; j < c; ++j)
      m.set(i, j, integer(0));
  return m;
}
RCP<const Basic> clean(const RCP<const Basic> &x) {
  return expand(simplify(x));
}
// Exact full-column-rank solve. Normal equations also handle MATLAB's
// overdetermined, consistent incidence systems without floating point.
DenseMatrix solve(const DenseMatrix &a, const DenseMatrix &b) {
  auto lhs = zeros(a.ncols(), a.ncols()), rhs = zeros(a.ncols(), b.ncols());
  if (a.nrows() == a.ncols()) {
    lhs = a;
    rhs = b;
  } else for (unsigned i = 0; i < a.ncols(); ++i) {
    for (unsigned j = 0; j < a.ncols(); ++j)
      for (unsigned k = 0; k < a.nrows(); ++k)
        lhs.set(i, j, add(lhs.get(i, j), mul(a.get(k, i), a.get(k, j))));
    for (unsigned j = 0; j < b.ncols(); ++j)
      for (unsigned k = 0; k < a.nrows(); ++k)
        rhs.set(i, j, add(rhs.get(i, j), mul(a.get(k, i), b.get(k, j))));
  }
  if (eq(*clean(det_berkowitz(lhs)), *integer(0)))
    throw std::runtime_error("singular QFA system");
  auto x = zeros(a.ncols(), b.ncols());
  fraction_free_gauss_jordan_solve(lhs, rhs, x);
  for (unsigned i = 0; i < x.nrows(); ++i)
    for (unsigned j = 0; j < x.ncols(); ++j)
      x.set(i, j, clean(x.get(i, j)));
  return x;
}
void weight(DenseMatrix &out, const DenseMatrix &v, unsigned row,
            const RCP<const Basic> &w) {
  for (unsigned i = 0; i < out.nrows(); ++i)
    for (unsigned j = 0; j < out.ncols(); ++j)
      out.set(i, j,
              add(out.get(i, j), mul(w, mul(v.get(row, i), v.get(row, j)))));
}
} // namespace
void complete_qfa(Topology &t, const ArchDef &arch) {
  const unsigned nc = t.N_caps, ns = t.N_sw, nn = arch.Acaps.nrows(),
                 no = t.m_ratios.nrows();
  for (unsigned i = 0; i < nc; ++i) {
    t.capacitances.push_back(symbol("C" + std::to_string(i + 1)));
    t.capacitor_esr.push_back(symbol("Resr" + std::to_string(i + 1)));
  }
  for (unsigned i = 0; i < ns; ++i)
    t.switch_resistances.push_back(symbol("Ron" + std::to_string(i + 1)));
  t.frequency = symbol("fsw");
  t.ZSSL = zeros(no, no);
  t.ZFSL = zeros(no, no);
  t.ZESR = zeros(no, no);
  for (auto &p : t.phase) {
    // Short supply. Capacitor nodal admittance enforces KCL and the
    // capacitor loop equations used by b_vec_multiphase.
    auto elements = zeros(p.graph.nrows(), p.graph.ncols() - 1),
         supply = zeros(p.graph.nrows(), 1);
    for (unsigned i = 0; i < p.graph.nrows(); ++i) {
      supply.set(i, 0, p.graph.get(i, 0));
      for (unsigned j = 1; j < p.graph.ncols(); ++j)
        elements.set(i, j - 1, p.graph.get(i, j));
    }
    auto sc = phase_conv(elements, supply);
    auto adm = zeros(sc.nrows(), sc.nrows()), loads = zeros(sc.nrows(), no);
    for (unsigned i = 0; i < sc.nrows(); ++i) {
      for (unsigned j = 0; j < sc.nrows(); ++j)
        for (unsigned k = 0; k < nc; ++k)
          adm.set(i, j,
                  add(adm.get(i, j),
                      mul(t.capacitances[k], mul(sc.get(i, k), sc.get(j, k)))));
      for (unsigned j = 0; j < no; ++j)
        loads.set(i, j, neg(sc.get(i, nc + j)));
    }
    auto potential = solve(adm, loads);
    p.b_vector = zeros(nc, no);
    p.r_vector = zeros(nc, no);
    for (unsigned k = 0; k < nc; ++k)
      for (unsigned j = 0; j < no; ++j) {
        RCP<const Basic> b = integer(0);
        for (unsigned i = 0; i < sc.nrows(); ++i)
          b = add(b, mul(t.capacitances[k],
                         mul(sc.get(i, k), potential.get(i, j))));
        p.b_vector.set(k, j, clean(b));
        p.r_vector.set(k, j, clean(sub(p.a_vector.get(k + 1, j), b)));
      }
    // MATLAB ar_builder: eliminate active-switch columns first, then
    // use only their pivot equations (not a least-squares KCL solve).
    auto q = zeros(nn, p.inc_on_conv_sw.ncols());
    for (unsigned i = 0; i < nn; ++i) {
      for (unsigned j = 0; j < p.n_on_sw; ++j)
        q.set(i, j, p.inc_on_sw.get(i, j));
      for (unsigned j = 0; j <= nc; ++j)
        q.set(i, p.n_on_sw + j, p.inc_on_conv_sw.get(i, j));
      for (unsigned j = 0; j < no; ++j)
        q.set(i, p.n_on_sw + nc + 1 + j,
              p.inc_on_conv_sw.get(i, nc + 1 + p.n_on_sw + j));
    }
    unsigned pivot_row = 0;
    for (unsigned col = 0; col < q.ncols() && pivot_row < nn; ++col) {
      unsigned pivot = pivot_row;
      while (pivot < nn && eq(*clean(q.get(pivot, col)), *integer(0)))
        ++pivot;
      if (pivot == nn) {
        if (col < p.n_on_sw)
          throw std::runtime_error("singular switch system");
        continue;
      }
      for (unsigned j = 0; j < q.ncols(); ++j) {
        auto v = q.get(pivot_row, j);
        q.set(pivot_row, j, q.get(pivot, j));
        q.set(pivot, j, v);
      }
      auto scale = q.get(pivot_row, col);
      for (unsigned j = 0; j < q.ncols(); ++j)
        q.set(pivot_row, j, div(q.get(pivot_row, j), scale));
      for (unsigned i = 0; i < nn; ++i)
        if (i != pivot_row) {
          auto factor = q.get(i, col);
          for (unsigned j = 0; j < q.ncols(); ++j)
            q.set(i, j,
                  clean(sub(q.get(i, j), mul(factor, q.get(pivot_row, j)))));
        }
      ++pivot_row;
    }
    auto sw = zeros(p.n_on_sw, no);
    for (unsigned i = 0; i < p.n_on_sw; ++i)
      for (unsigned j = 0; j < no; ++j) {
        auto v = q.get(i, p.n_on_sw + nc + 1 + j);
        for (unsigned k = 0; k <= nc; ++k)
          v = add(v, mul(q.get(i, p.n_on_sw + k), p.a_vector.get(k, j)));
        sw.set(i, j, clean(neg(v)));
      }
    p.ar_vector = zeros(nc + p.n_on_sw, no);
    for (unsigned i = 0; i < p.ar_vector.nrows(); ++i)
      for (unsigned j = 0; j < no; ++j)
        p.ar_vector.set(i, j,
                        i < nc ? p.a_vector.get(i + 1, j) : sw.get(i - nc, j));
    for (unsigned i = 0; i < nc; ++i) {
      weight(t.ZSSL, p.r_vector, i,
             div(integer(1),
                 mul(integer(2), mul(t.frequency, t.capacitances[i]))));
      weight(t.ZESR, p.ar_vector, i,
             div(t.capacitor_esr[i], p.duty.get_basic()));
    }
    for (unsigned i = 0; i < p.n_on_sw; ++i)
      weight(t.ZFSL, p.ar_vector, nc + i,
             div(t.switch_resistances[p.sw_idxs[i] - 1], p.duty.get_basic()));
  }
  t.ZSCC = zeros(no, no);
  for (unsigned i = 0; i < no; ++i)
    for (unsigned j = 0; j < no; ++j) {
      t.ZSSL.set(i, j, clean(t.ZSSL.get(i, j)));
      t.ZESR.set(i, j, clean(t.ZESR.get(i, j)));
      t.ZFSL.set(i, j, clean(add(t.ZFSL.get(i, j), t.ZESR.get(i, j))));
      t.ZSCC.set(i, j,
                 sqrt(add(pow(t.ZSSL.get(i, j), integer(2)),
                          pow(t.ZFSL.get(i, j), integer(2)))));
    }
  unsigned rows = 2 * (nc + 1) + t.phase[0].n_on_sw + t.phase[1].n_on_sw;
  auto a = zeros(rows, nc + 2 * nn), rhs = zeros(rows, 1);
  unsigned row = 0;
  for (unsigned ph = 0; ph < 2; ++ph) {
    a.set(row, nc + ph * nn, integer(1));
    rhs.set(row++, 0, integer(1));
    for (unsigned c = 0; c < nc; ++c, ++row) {
      a.set(row, c, integer(-1));
      for (unsigned n = 0; n < nn; ++n)
        a.set(row, nc + ph * nn + n, arch.Acaps.get(n, c));
    }
    for (unsigned c = 0; c < t.phase[ph].n_on_sw; ++c, ++row)
      for (unsigned n = 0; n < nn; ++n)
        a.set(row, nc + ph * nn + n, t.phase[ph].inc_on_sw.get(n, c));
  }
  auto volt = solve(a, rhs);
  t.vc = zeros(nc, 1);
  t.vr = zeros(2, ns);
  t.is = zeros(ns, no);
  for (unsigned c = 0; c < nc; ++c)
    t.vc.set(c, 0, volt.get(c, 0));
  for (unsigned ph = 0; ph < 2; ++ph) {
    for (unsigned c = 0; c < ns; ++c)
      for (unsigned n = 0; n < nn; ++n)
        t.vr.set(ph, c,
                 add(t.vr.get(ph, c),
                     mul(arch.Asw.get(n, c), volt.get(nc + ph * nn + n, 0))));
    for (unsigned c = 0; c < t.phase[ph].n_on_sw; ++c)
      for (unsigned j = 0; j < no; ++j)
        t.is.set(t.phase[ph].sw_idxs[c] - 1, j,
                 t.phase[ph].ar_vector.get(nc + c, j));
  }
}
