#include "solve_charge_vectors.h"
#include <symengine/integer.h>
#include <symengine/symbol.h>
#include <symengine/visitor.h>
#include <stdexcept>

using SymEngine::DenseMatrix;
using SymEngine::RCP;
using SymEngine::Basic;

static DenseMatrix zero_matrix(unsigned rows, unsigned cols) {
    DenseMatrix out(rows, cols);
    const auto zero = SymEngine::integer(0);
    for (unsigned r = 0; r < rows; ++r)
        for (unsigned c = 0; c < cols; ++c)
            out.set(r, c, zero);
    return out;
}

ChargeSolution solve_charge_vectors(const std::vector<DenseMatrix> &cutsets,
                                    unsigned n_caps, const DenseMatrix &duty,
                                    const std::vector<RCP<const SymEngine::Symbol>> &symbols) {
    if (cutsets.size() != 2) throw std::invalid_argument("solve_charge_vectors: exactly two phases required");
    if (!n_caps || duty.nrows() != 1 || duty.ncols() != 2)
        throw std::invalid_argument("solve_charge_vectors: invalid capacitor count or duty shape");
    for (unsigned c = 0; c < duty.ncols(); ++c)
        if (!duty.get(0, c).get())
            throw std::invalid_argument("solve_charge_vectors: null duty entry");
    const auto duty_symbols = SymEngine::free_symbols(*duty.get(0, 0));
    if (symbols.size() != duty_symbols.size())
        throw std::invalid_argument("solve_charge_vectors: symbols do not match duty");
    for (const auto &s : symbols)
        if (duty_symbols.count(s) == 0)
            throw std::invalid_argument("solve_charge_vectors: symbols do not match duty");

    if (cutsets[0].ncols() <= n_caps + 1)
        throw std::invalid_argument("solve_charge_vectors: output matrix empty");
    const unsigned outputs = cutsets[0].ncols() - (n_caps + 1);
    if (!SymEngine::eq(*duty.get(0, 1), *SymEngine::sub(SymEngine::integer(1), duty.get(0, 0))))
        throw std::invalid_argument("solve_charge_vectors: duty phases must sum to one");
    for (const auto &s : symbols)
        if (s.is_null()) throw std::invalid_argument("solve_charge_vectors: null symbol metadata");
    for (const auto &q : cutsets) {
        if (!q.nrows() || q.ncols() != n_caps + 1 + outputs)
            throw std::invalid_argument("solve_charge_vectors: inconsistent cutset dimensions");
        for (unsigned r = 0; r < q.nrows(); ++r)
            for (unsigned c = 0; c < q.ncols(); ++c)
                if (!q.get(r, c).get()) throw std::invalid_argument("solve_charge_vectors: null cutset entry");
    }
    const unsigned balance = n_caps;
    const unsigned system_size = 2 * (n_caps + 1);
    const unsigned cutset_rows = cutsets[0].nrows() + cutsets[1].nrows();
    if (cutset_rows < system_size - balance || cutset_rows > system_size)
        throw std::invalid_argument("solve_charge_vectors: cutset rows do not fill system");

    DenseMatrix qx = zero_matrix(system_size, system_size);
    DenseMatrix qo = zero_matrix(system_size, outputs);
    unsigned row = 0;
    std::vector<unsigned> phase_starts;
    for (unsigned p = 0; p < 2; ++p) {
        phase_starts.push_back(row);
        const DenseMatrix &q = cutsets[p];
        for (unsigned r = 0; r < q.nrows(); ++r, ++row) {
            qx.set(row, p * (n_caps + 1), SymEngine::neg(q.get(r, 0)));
            for (unsigned c = 1; c <= n_caps; ++c) qx.set(row, p * (n_caps + 1) + c, q.get(r, c));
            for (unsigned c = 0; c < outputs; ++c) qo.set(row, c, q.get(r, n_caps + 1 + c));
        }
    }
    for (unsigned c = 0; c < n_caps; ++c) {
        const unsigned r = system_size - balance + c;
        qx.set(r, c + 1, SymEngine::integer(1));
        qx.set(r, n_caps + 1 + c + 1, SymEngine::integer(1));
    }

    DenseMatrix rhs = zero_matrix(system_size, outputs);
    DenseMatrix ax = zero_matrix(system_size, outputs);
    for (unsigned r = 0; r < system_size; ++r)
        for (unsigned c = 0; c < outputs; ++c)
            rhs.set(r, c, SymEngine::neg(qo.get(r, c)));
    try {
        if (SymEngine::eq(*SymEngine::det_berkowitz(qx), *SymEngine::integer(0)))
            throw std::runtime_error("singular");
        SymEngine::fraction_free_LU_solve(qx, rhs, ax);
    } catch (...) {
        throw std::runtime_error("solve_charge_vectors: singular charge system");
    }

    ChargeSolution result{qx, qo, {}, zero_matrix(outputs, 1)};
    for (unsigned p = 0; p < 2; ++p) {
        DenseMatrix phase = zero_matrix(n_caps + 1, outputs);
        for (unsigned r = 0; r <= n_caps; ++r)
            for (unsigned c = 0; c < outputs; ++c) phase.set(r, c, ax.get(phase_starts[p] + r, c));
        result.a.push_back(phase);
        for (unsigned c = 0; c < outputs; ++c)
            result.m.set(c, 0, SymEngine::add(result.m.get(c, 0), phase.get(0, c)));
    }
    return result;
}
