#ifndef PHASE_CONV_H
#define PHASE_CONV_H

#include <symengine/matrix.h>

// Contract switches in MATLAB order. Aelem columns are [supply,caps,loads,off].
SymEngine::DenseMatrix phase_conv(const SymEngine::DenseMatrix &Aelem,
                                  const SymEngine::DenseMatrix &Asw);

#endif
