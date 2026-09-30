#include "dickson_matrix.h"
#include <vector>

void dickson_matrix(int n_stages, bool in_cap,
                    std::vector<std::vector<int>> &A_caps,
                    std::vector<std::vector<int>> &A_sw1,
                    std::vector<std::vector<int>> &A_sw2) {

    // Initialize the capacitor incidence matrix
    A_caps.assign(n_stages + 2, std::vector<int>(n_stages, 0));

    int n_sw1 = (n_stages > 2) ? 4 + (n_stages - 3) / 2 : 2;
    int n_sw2 = n_sw1 - (n_stages % 2);

    A_sw1.assign(n_stages + 3, std::vector<int>(n_sw1, 0));
    A_sw2.assign(n_stages + 3, std::vector<int>(n_sw2, 0));

    int j = 0;
    for (int i = 0; i < n_stages; ++i) {
        if (n_stages - (i + 1) < 3) {
            A_caps[i][i] = 1;
            if (i < n_stages - 1) {
                A_caps[n_stages + 1 - j][i] = -1;
                j++;
            }
        } else {
            A_caps[i][i] = 1;
            A_caps[i + 2][i] = -1;
        }
    }

    // Add the source line
    A_caps.insert(A_caps.begin(), std::vector<int>(n_stages, 0));

    for (int i = 0; i < n_sw2; ++i) {
        A_sw1[2 * i][i] = 1;
        if (i < n_sw2 - 1) {
            A_sw1[2 * i + 1][i] = -1;
        }
    }

    // Phase 2 Matrix
    for (size_t i = 1; i < A_sw2.size(); ++i) {
        for (size_t j = 0; j < A_sw2[0].size(); ++j) {
            A_sw2[i][j] = A_sw1[i - 1][j];
        }
    }

    if (n_stages > 2) {
        auto &special = (n_stages % 2 == 1) ? A_sw1 : A_sw2;
        const size_t row = special.size() - 3;
        const size_t col = special[0].size() - 2;
        for (size_t r = row; r < special.size(); ++r)
            for (size_t c = col; c < special[0].size(); ++c)
                special[r][c] = 0;
        special[row][col + 1] = 1;
        special[row + 1][col] = 1;
        special[row + 2][col + 1] = -1;
    } else {
        A_sw1.pop_back();
        A_sw2.pop_back();
        A_sw1.back().back() = -1;
        A_caps.erase(A_caps.end() - 2);
    }

    if (in_cap) {
        for (auto &row : A_caps) {
            row.insert(row.begin(), 0);
        }
        A_caps[0][0] = 1;
    }
}
