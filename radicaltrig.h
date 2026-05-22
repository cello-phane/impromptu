#ifndef RADICALTRIG_H
#define RADICALTRIG_H

#include <stdint.h>
#include <string.h> // memcpy for bitcast
#include <math.h> //for sqrt (or include SDL.h for SDL_sqrtf)

void rau_sincos(float in_rad, float *sin_out, float *cos_out);

// Convenience wrappers matching SDL signature
float rau_sinf(float x);
float rau_cosf(float x);
float rau_tanf(float x);

#endif
