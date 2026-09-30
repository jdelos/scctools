#include "dickson_arch.h"
#include "dickson_matrix.h"

#include <stdexcept>
#include <vector>

#include <symengine/integer.h>

ArchDef dickson_arch(int n_caps) {
    if (n_caps < 2) {
        throw std::invalid_argument("dickson_arch requires n_caps >= 2");
    }

    std::vector<std::vector<int>> caps, sw1, sw2;
    dickson_matrix(n_caps, false, caps, sw1, sw2);
    const unsigned n_switches = static_cast<unsigned>(sw1[0].size() + sw2[0].size());
    ArchDef arch{SymEngine::DenseMatrix(caps.size(), caps[0].size()),
                 SymEngine::DenseMatrix(sw1.size(), n_switches),
                 SymEngine::DenseMatrix(2, n_switches)};
    const auto zero = SymEngine::integer(0);
    for (unsigned r = 0; r < arch.Acaps.nrows(); ++r)
        for (unsigned c = 0; c < arch.Acaps.ncols(); ++c)
            arch.Acaps.set(r, c, zero);
    for (unsigned r = 0; r < arch.Asw.nrows(); ++r)
        for (unsigned c = 0; c < arch.Asw.ncols(); ++c)
            arch.Asw.set(r, c, zero);
    for (unsigned r = 0; r < arch.Asw_act.nrows(); ++r)
        for (unsigned c = 0; c < arch.Asw_act.ncols(); ++c)
            arch.Asw_act.set(r, c, zero);
    for (unsigned r = 0; r < caps.size(); ++r)
        for (unsigned c = 0; c < caps[r].size(); ++c)
            arch.Acaps.set(r, c, SymEngine::integer(caps[r][c]));
    for (unsigned r = 0; r < sw1.size(); ++r) {
        for (unsigned c = 0; c < sw1[r].size(); ++c) {
            arch.Asw.set(r, 2 * c, SymEngine::integer(sw1[r][c]));
            arch.Asw_act.set(0, 2 * c, SymEngine::integer(1));
        }
        for (unsigned c = 0; c < sw2[r].size(); ++c) {
            arch.Asw.set(r, 2 * c + 1, SymEngine::integer(sw2[r][c]));
            arch.Asw_act.set(1, 2 * c + 1, SymEngine::integer(1));
        }
    }
    return arch;
}
