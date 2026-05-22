#include "radicaltrig.h"

// Graph of sin/cos/tan and inverse functions - https://www.desmos.com/calculator/gkellct2v2

// Warp polynomial — maps t∈[0,1] → w∈[0,1], symmetric about 0.5
// ~5 digits accuracy
static const float C1=0.78539732f, C2=0.64607089f, C3=0.63401085f;
static const float C4=0.68518412f, C5=0.32482231f, C6=1.52006419f;

static float rwarp(float t) {
    float v  = t - 0.5f;
    float v2 = v * v;
    float p  = C1 + v2*(C2 + v2*(C3 + v2*(C4 + v2*(C5 + v2*C6))));
    return 0.5f + v * p;
}

// Radical Angle Unit - sincos, accurate to ~ 5 Digits
void rau_sincos(float in_rad, float *sin_out, float *cos_out) {
    // range reduction: fold into [0,1) RAU turns
    // radian to radical conversion - radian*2/π because 1 RAU == π/2
    //float rau = in_rad * (2.0f / 3.14159265358979323846f);
    float rau = in_rad; // input can be a radical angle unit if desired
    // shift to positive before floor
    float rau_pos = rau + 4096.0f;        // bias to ensure positive
    // Since 4096 = 2¹² and 4096 & 3 = 0, adding it shifts the integer
    // but leaves qi = qi_full & 3 unchanged relative to what it would be without the bias
    int   quadrant_index_full = (int)rau_pos;  // integral is an index to map result
    float frac = rau_pos - (float)quadrant_index_full; // fractional ∈ [0,1)
    int   quadrant_index = quadrant_index_full & 3; // quadrant 0 to 3

    // warp t into w so linear --> arc length (uniformly mapped)
    float w = rwarp(frac);

    // odd-quadrant reversal: Q1,Q3 → w = 1-w
    if (quadrant_index & 1) w = 1.0f - w;

    // Diagonal to unit circle
    float omw = 1.0f - w;
    float D   = omw*omw + w*w;    // = 1 - 2w(1-w)
    float inv = 1.0f / sqrtf(D);  // Trig formulae:
    float cs  = omw * inv;        // (1-w) ÷ sqrt(1 + 2×w + 2×w^2)
    float sn  = w   * inv;        //  w ÷ sqrt(1 + 2×w + 2×w^2)

    // sign bits from quadrant
    // cos negative in Q1,Q2 (qi=1,2): csign bit set when (qi+1)>>1 & 1
    // sin negative in Q2,Q3 (qi=2,3): ssign bit set when  qi>>1    & 1
    uint32_t csign = (uint32_t)(((quadrant_index+1)>>1) & 1) << 31;
    uint32_t ssign = (uint32_t)( (quadrant_index>>1)    & 1) << 31;

    uint32_t cs_bits, sn_bits;
    memcpy(&cs_bits, &cs, 4);
    memcpy(&sn_bits, &sn, 4);
    cs_bits ^= csign;
    sn_bits ^= ssign;
    memcpy(cos_out, &cs_bits, 4);
    memcpy(sin_out, &sn_bits, 4);
}

// Convenience wrappers
float rau_sinf(float x) {
    float s, c;
    rau_sincos(x, &s, &c);
    return s;
}

float rau_cosf(float x) {
    float s, c; rau_sincos(x, &s, &c);
    return c;
}

float rau_tanf(float x) {
    float s, c; rau_sincos(x, &s, &c);
    // this version — acceptable for games/graphics where
    // a missed spike is better than a NaN if user knows it's determined?
    return (fabsf(c) < 1e-6f) ? 0.0f : s / c;
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
