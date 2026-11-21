// Public functions with no internal counterpart: the version and accessors for the opaque types.
#include "expression.h"

// PS_VERSION of the built library (to compare with the header's at run time)
const char* ps_version(void) { return PS_VERSION; }

// e as a polynomial (borrowed, valid while e lives), or NULL when e is not a polynomial
const ps_poly* ps_expr_poly(const ps_expr* e) { return e && e->type == EXPR_POLY ? e->poly : NULL; }
