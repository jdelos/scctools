#ifndef QFA_BUNDLE_H
#define QFA_BUNDLE_H
#include "dickson_hybrid_topology.h"
// Complete graph QFA, following generic_switched_capacitor_class.m and
// SCC_Phase.m (Julia Delos). ZSCC is the elementwise analytical RSS
// approximation.
void complete_qfa(Topology &top, const ArchDef &arch);
#endif
