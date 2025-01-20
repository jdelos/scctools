#ifndef DICKSON_MATRIX_H
#define DICKSON_MATRIX_H

#include <vector>

// Function to generate the incidence matrices for the Dickson Ladder converter
void dickson_matrix(int n_stages, bool in_cap,
                    std::vector<std::vector<int>> &A_caps,
                    std::vector<std::vector<int>> &A_sw1,
                    std::vector<std::vector<int>> &A_sw2);

#endif // DICKSON_MATRIX_H

