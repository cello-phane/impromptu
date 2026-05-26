#ifndef RADICALTRIG_H
#define RADICALTRIG_H

#include <SDL.h>
//#include <stdint.h>
//#include <string.h> // memcpy for bitcast
//#include <math.h> //for sqrt (or include SDL.h for SDL_sqrtf)

// prototypes
void rau_sincos(float input_t, float *sin_out, float *cos_out);
float rau_sinf(float x);
float rau_cosf(float x);
float rau_tanf(float x);

#endif
