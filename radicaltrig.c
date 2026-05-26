#include "radicaltrig.h"
#include "SDL_stdinc.h"
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

static unsigned int float_to_bits(float x)
{
    unsigned int u;
    SDL_memcpy(&u, &x, sizeof u);
    return u;
}

static float bits_to_float(unsigned int u)
{
    float x;
    SDL_memcpy(&x, &u, sizeof x);
    return x;
}

// Radical Angle Unit - sincos, accurate to ~ 5 Digits
void rau_sincos(float input_t, float *sin_out, float *cos_out) {
    // if radian is input type and NON_UNIFORM_VEL is 0:
    // input_t = input_t * (2/π);
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
    unsigned int csign = (unsigned int)(((quadrant_index+1)>>1) & 1) << 31;
    unsigned int ssign = (unsigned int)( (quadrant_index>>1)    & 1) << 31;

    unsigned int cs_bits = float_to_bits(cs) ^ csign;
    unsigned int sn_bits = float_to_bits(sn) ^ ssign;

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
    // if radian is input type and NON_UNIFORM_VEL is 0:
    // x = x * (2/π);
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
