#include "dickson_arch.h"
#include "dickson_matrix.h"
#include <cassert>
#include <stdexcept>

int main() {
    std::vector<std::vector<int>> caps, sw1, sw2;
    dickson_matrix(2, false, caps, sw1, sw2);
    assert(caps.size() == 3 && caps.front().size() == 2);
    assert(sw1.size() == 5 && sw1.front().size() == 2);
    assert(sw2.size() == 5 && sw2.front().size() == 2);
    auto arch = dickson_arch(2);
    assert(arch.A.nrows() == 3 && arch.A.ncols() == 2);
    assert(arch.m.nrows() == 1 && arch.m.ncols() == 1);
    bool rejected = false;
    try { dickson_matrix(1, false, caps, sw1, sw2); }
    catch (const std::invalid_argument &) { rejected = true; }
    assert(rejected);
}
