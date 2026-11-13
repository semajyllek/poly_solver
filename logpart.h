#ifndef LOGPART_H
#define LOGPART_H

#include "expression.h"

// Logarithmic part of ∫ a/d by Lazard–Rioboo–Trager, for square-free d with
// deg a < deg d and deg d >= 2. Returns ln/atan terms and, when a log coefficient is an
// algebraic number of degree >= 3, an exact RootSum. Returns NULL only when a == 0.
Expr* lrt_log_part(const Polynomial* a, const Polynomial* d);

#endif
