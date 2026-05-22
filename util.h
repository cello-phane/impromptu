#ifndef UTIL_H
#define UTIL_H

#include <string.h>

#include <SDL.h>

#define MAX(a, b) ((a > b) ? a : b)
#define MIN(a, b) ((a < b) ? a : b)

// math.h gates M_PI behind feature-test macros under strict C99.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

Uint32 rand32(Uint64 *rng);

#endif