#include "radicaltrig.h"
#include "SDL_stdinc.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

// Graph of sin/cos/tan and inverse functions - https://www.desmos.com/calculator/gkellct2v2

// NON_UNIFORM_VEL 0 — arc-uniform: warp polynomial applied, constant arc speed (π/2 per unit parameter)
// NON_UNIFORM_VEL 1 — raw linear diagonal projection, no warp, fast, non-uniform in arc
#define NON_UNIFORM_VEL 0

// ── Helper Functions ────────────────────────────────────────
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

// ── Forward Trigonometric Functions ────────────────────────────────────────
float rau_warpf(float t) {
    static const float C[6] = {
        0.78539732f,   /* ≈ π/4  — Remez-confirmed leading coefficient */
        0.64607089f,
        0.63401085f,
        0.68518412f,
        0.32482231f,
        1.52006419f
    };
    float v  = t - 0.5f;    /* centre on 0 to exploit odd symmetry */
    float v2 = v * v;
    float p  = C[5];        /* Nested sum of products in v² (even powers only) */
    for (int i = 4; i >= 0; --i) p = v2 * p + C[i];
    return 0.5f + v * p;
}

void rau_sincosf(float input_t, float *sin_out, float *cos_out) {
    /* Range reduction: fold into [0,4), extract integer quadrant */
    float rau_pos = mod4(input_t);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;  /* fractional part ∈ [0,1) */
    int   qi      = qi_full & 3;               /* quadrant: 0,1,2,3 */

    /* Arc-length correction (compile-time switch) */
#if NON_UNIFORM_VEL
    float w = frac;                /* raw diagonal: non-uniform arc speed */
#else
    float w = rau_warpf(frac);     /* warped: uniform arc speed = π/2 per RAU */
#endif

    /* Odd-quadrant reversal: w increases [0→1] in all quadrants after this */
    if (qi & 1) w = 1.0f - w;

    /* Diagonal point to unit circle — single sqrt, both outputs */
    float omw = 1.0f - w;
    float D   = omw*omw + w*w;         /* squared distance from origin to diagonal point */
    float inv = 1.0f / SDL_sqrtf(D);   /* normalisation factor */
    float cs  = omw * inv;             /* cos = (1-w)/sqrt(D) */
    float sn  = w   * inv;             /* sin =    w /sqrt(D) */

    Uint32 csign = (Uint32)(((qi+1)>>1) & 1) << 31;
    Uint32 ssign = (Uint32)( (qi>>1)    & 1) << 31;

    Uint32 cs_bits = float_to_bits(cs) ^ csign;
    Uint32 sn_bits = float_to_bits(sn) ^ ssign;

    *cos_out = bits_to_float(cs_bits);
    *sin_out = bits_to_float(sn_bits);
}

float rau_sinf(float x) {
    float s, c;
    rau_sincosf(x, &s, &c);
    return s;
}

float rau_cosf(float x) {
    float s, c;
    rau_sincosf(x, &s, &c);
    return c;
}

float rau_tanf(float x) {
    float rau_pos = mod4(x);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;
    int   qi      = qi_full & 3;

#if NON_UNIFORM_VEL
    float w = frac;
#else
    float w = rau_warpf(frac);
#endif

    if (qi & 1) w = 1.0f - w;

    float denom = 1.0f - w;               /* = cos (without normalisation) */
    if (denom < 1e-6f) denom = 1e-6f;    /* pole clamp — see note above */

    float t = w / denom;                  /* = tan (sqrt cancels) */
    if (qi & 1) t = -t;                   /* sign: tan negative in Q1,Q3 */
    return t;
}

void rau_sincos_mf(float input_t, float M, float *sin_out, float *cos_out) {
    float rau_pos = mod4(input_t);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;
    int   qi      = qi_full & 3;

    /* Lerp between raw (M=0) and warped (M=1) */
    float w_raw  = frac;
    float w_warp = rau_warpf(frac);
    float w      = w_raw + M * (w_warp - w_raw);

    if (qi & 1) w = 1.0f - w;

    float omw = 1.0f - w;
    float D   = omw*omw + w*w;
    float inv = 1.0f / SDL_sqrtf(D);
    float cs  = omw * inv;
    float sn  = w   * inv;

    Uint32 csign = (Uint32)(((qi+1)>>1) & 1) << 31;
    Uint32 ssign = (Uint32)( (qi>>1)    & 1) << 31;

    *cos_out = bits_to_float(float_to_bits(cs) ^ csign);
    *sin_out = bits_to_float(float_to_bits(sn) ^ ssign);
}

float rau_arctan_adjf(float x) {
    float x2 = x * x;
    return x * (
        1.000087f
        + x2 * (-0.33288950512027f
        + x2 * ( 0.19383271707398f
        + x2 * (-0.11735031947869f
        + x2 * ( 0.05368137843104f
        + x2 * (-0.01213232131734f)))))
    );
}

float rau_r_arctanf(float ry, float rx, int *err) {
    if (err) *err = 0;
    if (ry == 0.0f && rx == 0.0f) {
        if (err) *err = 1;
        return NAN;
    }
    if (rx == 0.0f) return 1.0f;         /* pure sin-axis: w = 1 */
    float t = SDL_fabsf(ry / rx);
    return t / (1.0f + t);               /* w = tan/(1+tan) — exact rational */
}

float rau_r_arcsinf(float s, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(s) || s < -1.0f || s > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float sy2   = s * s;
    float denom = 2.0f * sy2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f) return 0.5f;   /* 45° singularity guard */
    float disc  = SDL_sqrtf(SDL_max(sy2 * (1.0f - sy2), 0.0f));
    return SDL_fabsf((sy2 - disc) / denom);
}

float rau_r_arccosf(float c, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(c) || c < -1.0f || c > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float cx2   = c * c;
    float denom = 2.0f * cx2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f) return 0.5f;   /* 45° singularity guard */
    float disc  = SDL_sqrtf(SDL_max(cx2 * (1.0f - cx2), 0.0f));
    return SDL_fabsf((cx2 - 1.0f + disc) / denom);
}

float rau_invpolyf(float w, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(w) || w < 0.0f || w > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    if (w >= 0.5f) {
        /* reflect: use complementary angle identity π/2 - arctan(1/u) = arctan(u) */
        float u = (1.0f - w) / w;
        float a = u * (
            1.000087f
            + (u*u) * (-0.33288950512027f
            + (u*u) * ( 0.19383271707398f
            + (u*u) * (-0.11735031947869f
            + (u*u) * ( 0.05368137843104f
            + (u*u) * (-0.01213232131734f)))))
        );
        return ((float)M_PI_2 - a) * (2.0f / (float)M_PI);
    } else {
        float u = w / (1.0f - w);
        return rau_arctan_adjf(u) * (2.0f / (float)M_PI);
    }
}

float rau_atan2f(float ry, float rx, int *err) {
    if (err) *err = 0;

    if (!rau_isfinitef(ry) || !rau_isfinitef(rx)) {
        if (err) *err = 1;
        return NAN;
    }
    if (ry == 0.0f && rx == 0.0f) {
        if (err) *err = 1;
        return NAN;
    }

    float w = rau_r_arctanf(ry, rx, err);
    if (err && *err) return NAN;

    float t = rau_invpolyf(w, err);
    if (err && *err) return NAN;

    /* Quadrant offset + fractional phase.
     * Q0 [0°- 90°]: t increases forward   (0 + t)
     * Q1 [90°-180°]: t increases reversed (1 + (1-t))
     * Q2 [180°-270°]: t increases forward (2 + t)
     * Q3 [270°-360°]: t increases reversed (3 + (1-t))
     */
    if (rx >= 0.0f && ry >= 0.0f) return 0.0f + t;
    if (rx <  0.0f && ry >= 0.0f) return 1.0f + (1.0f - t);
    if (rx <  0.0f && ry <  0.0f) return 2.0f + t;
    return                                 3.0f + (1.0f - t);
}
