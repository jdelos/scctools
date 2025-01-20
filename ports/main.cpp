#include "dickson_hybrid_topology.h"
#include <iostream>

int main() {
    Topology topology = dickson_hybrid_topology(5, 0.5, true, false);
    std::cout << "Topology calculated for " << topology.N_caps << " capacitors." << std::endl;
    return 0;
}
