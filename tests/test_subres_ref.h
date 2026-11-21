#ifndef TEST_SUBRES_REF_H
#define TEST_SUBRES_REF_H

#include "subres.h"

// S_j(t, x) by exact Sylvester determinants and interpolation (independent test oracle)
PolyTX ref_subresultant(const RTPair* p, int j);

#endif
