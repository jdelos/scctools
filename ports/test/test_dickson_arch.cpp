#include "dickson_arch.h"
#include <symengine/integer.h>
#include <symengine/basic.h>
#include <cassert>
#include <stdexcept>

static void check(int n, const int *caps, unsigned cap_rows, unsigned cap_cols,
                  const int *sw, unsigned sw_rows, unsigned sw_cols) {
    const ArchDef a = dickson_arch(n);
    assert(a.Acaps.nrows() == cap_rows && a.Acaps.ncols() == cap_cols);
    assert(a.Asw.nrows() == sw_rows && a.Asw.ncols() == sw_cols);
    assert(a.Asw_act.nrows() == 2 && a.Asw_act.ncols() == sw_cols);
    for (unsigned r = 0; r < cap_rows; ++r)
        for (unsigned c = 0; c < cap_cols; ++c)
            assert(SymEngine::eq(*a.Acaps.get(r, c), *SymEngine::integer(caps[r * cap_cols + c])));
    for (unsigned r = 0; r < sw_rows; ++r)
        for (unsigned c = 0; c < sw_cols; ++c)
            assert(SymEngine::eq(*a.Asw.get(r, c), *SymEngine::integer(sw[r * sw_cols + c])));
    for (unsigned c = 0; c < sw_cols; ++c)
        for (unsigned r = 0; r < 2; ++r)
            assert(SymEngine::eq(*a.Asw_act.get(r, c), *SymEngine::integer(r == c % 2)));
}

int main() {
    const int caps2[] = {0,0,1,0,0,1,-1,0};
    const int sw2[] = {1,0,0,0,-1,1,0,0,0,-1,1,0,0,0,-1,1};
    check(2, caps2, 4, 2, sw2, 4, 4);
    const int caps3[] = {0,0,0,1,0,0,0,1,0,0,0,1,0,-1,0,-1,0,0};
    const int sw3[] = {1,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,-1,1,0,0,0,0,0,0,-1,1,0,0,1,0,0,0,-1,1,0,0,0,0,0,0,0,1,-1};
    check(3, caps3, 6, 3, sw3, 6, 7);
    const ArchDef a4 = dickson_arch(4);
    assert(a4.Acaps.nrows() == 7 && a4.Acaps.ncols() == 4);
    assert(a4.Asw.nrows() == 7 && a4.Asw.ncols() == 8);
}
