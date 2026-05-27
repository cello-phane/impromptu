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

float rau_arctan_adj(float x);
float rau_r_arctan(float ry, float rx, int *err);
float rau_r_arcsin(float s, int *err);
float rau_r_arccos(float c, int *err);
float rau_invpoly(float w, int *err);
float rau_full_phi(float ry, float rx, int *err);

#endif
