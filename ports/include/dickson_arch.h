#ifndef DICKSON_ARCH_H
#define DICKSON_ARCH_H

#include <symengine/matrix.h>

struct ArchDef {
    SymEngine::DenseMatrix Acaps;
    SymEngine::DenseMatrix Asw;
    SymEngine::DenseMatrix Asw_act;
};

ArchDef dickson_arch(int n_caps);

#endif

