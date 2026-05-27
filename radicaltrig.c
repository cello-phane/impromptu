#include "radicaltrig.h"
#include "SDL_stdinc.h"
#include <math.h>

// Graph of sin/cos/tan and inverse functions - https://www.desmos.com/calculator/gkellct2v2

// NON_UNIFORM_VEL 0 — arc-uniform: warp polynomial applied,
//                     constant arc speed (π/2 per unit parameter)
//                     correct for FFT twiddles, motor control, camera pans
//
// NON_UNIFORM_VEL 1 — raw linear diagonal projection, no warp,
//                     faster (one fewer polynomial evaluation + sqrt),
//                     non-uniform arc speed (peaks 2× at 45°),
//                     acceptable for rendering, easing, visual effects
#define NON_UNIFORM_VEL 0

// Warp polynomial — maps t∈[0,1] → w∈[0,1], symmetric about 0.5 [~5 digits accuracy]
static const float C1=0.78539732f, C2=0.64607089f, C3=0.63401085f;
static const float C4=0.68518412f, C5=0.32482231f, C6=1.52006419f;
// Approximates sin(θ)/(sin(θ)+cos(θ)) for diagonal to arc length stretching
static float rwarp(float t) {
    float v  = t - 0.5f;
    float v2 = v * v;
    float p  = C1 + v2*(C2 + v2*(C3 + v2*(C4 + v2*(C5 + v2*C6))));
    return 0.5f + v * p;
}

// Helpers
static float mod4(float a) {
    float r = SDL_fmodf(a, 4.0f);
    if (r < 0.0f) r += 4.0f;
    return r;
}

static Uint32 float_to_bits(float x) {
    Uint32 u; SDL_memcpy(&u, &x, sizeof u); return u;
}

static float bits_to_float(Uint32 u) {
    float x; SDL_memcpy(&x, &u, sizeof x); return x;
}

static inline int rau_isfinitef(float x) {
    return isfinite(x);
}

static inline int rau_in_unit_range(float x) {
    return x >= -1.0f && x <= 1.0f;
}

static inline float rau_clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

static inline float rau_sanitize_unit(float x) {
    if (x < -1.0f) return -1.0f;
    if (x >  1.0f) return  1.0f;
    return x;
}

// Radical Angle Unit - sincos, accurate to ~ 5 Digits
void rau_sincos(float input_t, float *sin_out, float *cos_out) {
    // Input units: RAU (0 to 4 = full revolution)
    // For radian input: input_t *= (2.0f / M_PI)  before calling

    /* --- Range reduction: fold into [0,1) and integer is a quadrant index --- */
    float rau_pos = mod4(input_t);
    int   quadrant_index_full = (int)rau_pos;
    float frac = rau_pos - (float)quadrant_index_full;
    int   quadrant_index = quadrant_index_full & 3; // quadrant 0 to 3

    // Warp into w making linear --> arc length
    #if NON_UNIFORM_VEL
    float w = frac;
    #else
    float w = rwarp(frac);
    #endif

    // Odd-quadrant reversal: Q1,Q3 → w = 1-w
    if (quadrant_index & 1) w = 1.0f - w;

    // Diagonal to unit circle
    float omw = 1.0f - w;
    float D   = omw*omw + w*w;    // = 1 - 2w(1-w)
    float inv = 1.0f / SDL_sqrtf(D);  // Trig formulae:
    float cs  = omw * inv;        // (1-w) ÷ sqrt(1 - 2w + 2w^2)
    float sn  = w   * inv;        //     w ÷ sqrt(1 - 2w + 2w^2)

    // Sign bits from quadrant:
    // cos negative in Q1,Q2 (qi=1,2): csign bit set when (qi+1)>>1 & 1
    // sin negative in Q2,Q3 (qi=2,3): ssign bit set when  qi>>1    & 1
    Uint32 csign = (Uint32)(((quadrant_index+1)>>1) & 1) << 31;
    Uint32 ssign = (Uint32)( (quadrant_index>>1)    & 1) << 31;

    Uint32 cs_bits = float_to_bits(cs) ^ csign;
    Uint32 sn_bits = float_to_bits(sn) ^ ssign;

    *cos_out = bits_to_float(cs_bits);
    *sin_out = bits_to_float(sn_bits);
}

// Convenience wrappers
float rau_sinf(float x) {
    float s, c;
    rau_sincos(x, &s, &c);
    return s;
}

float rau_cosf(float x) {
    float s, c;
    rau_sincos(x, &s, &c);
    return c;
}

// Standalone function
float rau_tanf(float x)
{
    // Input units: RAU (0 to 4 = full revolution)
    // For radian input: x *= (2.0f / M_PI)  before calling
    float rau_pos = mod4(x);
    int   quadrant_index_full = (int)rau_pos;
    float frac = rau_pos - (float)quadrant_index_full;
    int   quadrant_index = quadrant_index_full & 3; // quadrant 0 to 3

    // warp into w making linear --> arc length
    #if NON_UNIFORM_VEL
    float w = frac;
    #else
    float w = rwarp(frac);
    #endif

    if (quadrant_index & 1)
        w = 1.0f - w;
    float denom = 1.0f - w;

    // optional clamp
    if (denom < 1e-6f)
        denom = 1e-6f;
    float t = w / denom;
    if (quadrant_index & 1)
        t = -t;

    return t;
}

// Hot path version with M to set velocity
void rau_sincos_m(float input_t, float M, float *sin_out, float *cos_out) {
    // M=0: NON_UNIFORM_VEL behaviour (linear diagonal)
    // M=1: arc-uniform behaviour (warp polynomial)
    // M in between: smooth interpolation — CVT blend
    float rau_pos = mod4(input_t);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;
    int   qi      = qi_full & 3;
    float w_raw   = frac;
    float w_warp  = rwarp(frac);
    float w       = w_raw + M * (w_warp - w_raw);  // lerp

    // Odd-quadrant reversal: Q1,Q3 → w = 1-w
    if (qi & 1) w = 1.0f - w;

    // Diagonal to unit circle
    float omw = 1.0f - w;
    float D   = omw*omw + w*w;    // = 1 - 2w(1-w)
    float inv = 1.0f / SDL_sqrtf(D);  // Trig formulae:
    float cs  = omw * inv;        // (1-w) ÷ sqrt(1 - 2w + 2w^2)
    float sn  = w   * inv;        //     w ÷ sqrt(1 - 2w + 2w^2)

    // Sign bits from quadrant:
    // cos negative in Q1,Q2 (qi=1,2): csign bit set when (qi+1)>>1 & 1
    // sin negative in Q2,Q3 (qi=2,3): ssign bit set when  qi>>1    & 1
    Uint32 csign = (Uint32)(((qi+1)>>1) & 1) << 31;
    Uint32 ssign = (Uint32)( (qi>>1)    & 1) << 31;

    Uint32 cs_bits = float_to_bits(cs) ^ csign;
    Uint32 sn_bits = float_to_bits(sn) ^ ssign;

    *cos_out = bits_to_float(cs_bits);
    *sin_out = bits_to_float(sn_bits);
}

float rau_arctan_adj(float x) {
    float x2 = x * x;
    return x * (
        1.000087f
        + x2 * (-0.33288950512027f
        + x2 * (0.19383271707398f
        + x2 * (-0.11735031947869f
        + x2 * (0.05368137843104f
        + x2 * (-0.01213232131734f)))))
    );
}

float rau_r_arctan(float ry, float rx, int *err) {
    if (err) *err = 0;
    if (ry == 0.0f && rx == 0.0f) {
        if (err) *err = 1;
        return NAN;
    }
    if (rx == 0.0f) return 1.0f;
    float t = SDL_fabsf(ry / rx);
    return t / (1.0f + t);
}

float rau_r_arcsin(float s, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(s) || s < -1.0f || s > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float sy2 = s * s;
    float denom = 2.0f * sy2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f) return 0.5f;
    float disc = SDL_sqrtf(SDL_max(sy2 * (1.0f - sy2), 0.0f));
    return SDL_fabsf((sy2 - disc) / denom);
}

float rau_r_arccos(float c, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(c) || c < -1.0f || c > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float cx2 = c * c;
    float denom = 2.0f * cx2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f) return 0.5f;
    float disc = SDL_sqrtf(SDL_max(cx2 * (1.0f - cx2), 0.0f));
    return SDL_fabsf((cx2 - 1.0f + disc) / denom);
}

float rau_invpoly(float w, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(w) || w < 0.0f || w > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    static const float AC[6] = {
        1.000087f, -0.33288950512027f, 0.19383271707398f,
        -0.11735031947869f, 0.05368137843104f, -0.01213232131734f
    };
    if (w >= 0.5f) {
        float u = (1.0f - w) / w;
        return (float)(M_PI_2 - (u * (
            AC[0] + (u*u) * (AC[1] + (u*u) * (AC[2] + (u*u) * (AC[3] + (u*u) * (AC[4] + (u*u) * AC[5]))))
        ))) * (2.0f / (float)M_PI);
    } else {
        float u = w / (1.0f - w);
        return rau_arctan_adj(u) * (2.0f / (float)M_PI);
    }
}

float rau_full_phi(float ry, float rx, int *err) {
    if (err) *err = 0;
    if (ry == 0.0f && rx == 0.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float w = rau_r_arctan(ry, rx, err);
    if (w==0.0) return NAN;
    float t = rau_invpoly(w, err);
    if (t==0.0) return NAN;

    if (rx >= 0.0f && ry >= 0.0f) return 0.0f + t;
    if (rx <  0.0f && ry >= 0.0f) return 1.0f + (1.0f - t);
    if (rx <  0.0f && ry <  0.0f) return 2.0f + t;
    return 3.0f + (1.0f - t);
}

/* LICENSE
 [Michael Ledesma - cello-phane@github]
 This is free and unencumbered software released into the public domain.

 Anyone is free to copy, modify, publish, use, compile, sell, or
 distribute this software, either in source code form or as a compiled
 binary, for any purpose, commercial or non-commercial, and by any
 means.

 In jurisdictions that recognize copyright laws, the author or authors
 of this software dedicate any and all copyright interest in the
 software to the public domain. We make this dedication for the benefit
 of the public at large and to the detriment of our heirs and
 successors. We intend this dedication to be an overt act of
 relinquishment in perpetuity of all present and future rights to this
 software under copyright law.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 OTHER DEALINGS IN THE SOFTWARE.

 For more information, please refer to <https://unlicense.org>

 */
