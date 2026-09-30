#include "dickson_matrix.h"
#include <cassert>
#include <vector>

int main() {
    std::vector<std::vector<int>> A_caps, A_sw1, A_sw2;
    dickson_matrix(3, true, A_caps, A_sw1, A_sw2);

    assert(A_caps.size() == 6);
    assert(A_caps[0].size() == 4);
    assert(A_caps[0][0] == 1);
    assert(A_caps[1][1] == 1);
    assert(A_sw1.size() == 6);
    assert(A_sw2.size() == 6);
}
