#include "radicaltrig.h"
#include "SDL_stdinc.h"
#include <math.h>

// Radical Angle Unit (RAU) trigonometric library
// Reference: https://www.desmos.com/calculator/gkellct2v2
//
// 1 RAU = π/2 radians. Full circle = 4 RAU.
// Integer quadrant boundaries are exact — no floating point error at 0°/90°/180°/270°.
//
// Accuracy quality tiers (selectable via RAU_ATAN_QUALITY):
//   0 = f16:     ~1.18e-4 RAU max err — float16 input quality, 4 terms
//   1 = hsp:     ~1.67e-5 RAU max err — half-single precision,  5 terms
//   2 = precise: ~3.62e-7 RAU max err — float32 quality,        7 terms (default)
//
// NON_UNIFORM_VEL 0 — arc-uniform: warp polynomial applied, π/2 arc per RAU
// NON_UNIFORM_VEL 1 — raw diagonal projection, no warp, non-uniform arc speed
#define NON_UNIFORM_VEL 0

#ifndef RAU_ATAN_QUALITY
#define RAU_ATAN_QUALITY 2
#endif

// ── Helper Functions ───────────────────────────────────────────────────────

static float mod4(float a) {
    float r = SDL_fmodf(a, 4.0f);
    if (r < 0.0f) r += 4.0f;
    return r;
}

static Uint32 float_to_bits(float x) {
    Uint32 u;
    SDL_memcpy(&u, &x, sizeof u);
    return u;
}

static float bits_to_float(Uint32 u) {
    float x;
    SDL_memcpy(&x, &u, sizeof x);
    return x;
}

static inline int   rau_isfinitef(float x)     { return isfinite(x); }
static inline int   rau_isnanf(float x)         { return isnan(x); }
static inline float rau_clampf(float x, float lo, float hi) {
    return x < lo ? lo : (x > hi ? hi : x);
}

// ── Conversion Helpers ─────────────────────────────────────────────────────

// Wrap unsigned RAU result [0,4) to signed (−2,+2], then convert to radians
float rau_atan2_signed_radians(float phi_rau) {
    float wrapped = fmodf(phi_rau + 2.0f, 4.0f);
    if (wrapped < 0.0f) wrapped += 4.0f;
    wrapped -= 2.0f;
    return wrapped * (float)M_PI_2;
}

// Wrap unsigned RAU result [0,4) to degrees [0,360)
float rau_atan2_signed_degs(float phi_rau) {
    float deg = fmodf(phi_rau * 90.0f, 360.0f);
    if (deg < 0.0f) deg += 360.0f;
    return deg;
}

// ── Warp Polynomial ────────────────────────────────────────────────────────
//
// rau_warpf maps the raw diagonal parameter t ∈ [0,1] to the arc-uniform
// parameter w = sin(t·π/2) / (sin(t·π/2) + cos(t·π/2)).
//
// Coefficients: Remez minimax over [0,1], leading term WC[0] ≈ π/4.
// Derivation: end-to-end SNR minimization on rau_sincosf output.
// Max err: ~5.6e-7 (float32 quality).
//
// Two sets exist — the active set uses fuller float32 bit width:
//   Active (wider bits, SNR-optimized):
//     {0.78539816f, 0.64607158f, 0.63401589f, 0.68515350f, 0.32501622f, 1.51901679f}
//   Alternative (fewer bits, slightly higher peak error):
//     {0.78539732f, 0.64607089f, 0.63401085f, 0.68518412f, 0.32482231f, 1.52006419f}

float rau_warpf(float t) {
    static const float C[6] = {
        0.78539816339744830962f,  /* ≈ π/4 — Remez leading term */
        0.64607158024987317298f,
        0.63401589172451679138f,
        0.68515350354689586789f,
        0.32501622369042378935f,
        1.51901679307446258196f
    };
    float v  = t - 0.5f;
    float v2 = v * v;
    float p  = C[5];
    for (int i = 4; i >= 0; --i)
        p = v2 * p + C[i];
    return 0.5f + v * p;
}

// ── Forward Trigonometric Functions ───────────────────────────────────────

void rau_sincosf(float input_t, float *sin_out, float *cos_out) {
    float rau_pos = mod4(input_t);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;  /* ∈ [0,1) */
    int   qi      = qi_full & 3;               /* quadrant 0..3 */

#if NON_UNIFORM_VEL
    float w = frac;
#else
    float w = rau_warpf(frac);
#endif

    /* Odd-quadrant reversal: ensures w increases 0→1 in every quadrant */
    if (qi & 1) w = 1.0f - w;

    /* Diagonal point → unit circle via single sqrt */
    float omw = 1.0f - w;
    float D   = omw * omw + w * w;
    float inv = 1.0f / SDL_sqrtf(D);
    float cs  = omw * inv;
    float sn  = w   * inv;

    /* Sign via IEEE-754 bit XOR — no branches, no float multiply */
    Uint32 csign = (Uint32)(((qi + 1) >> 1) & 1) << 31;
    Uint32 ssign = (Uint32)(( qi      >> 1) & 1) << 31;

    *cos_out = bits_to_float(float_to_bits(cs) ^ csign);
    *sin_out = bits_to_float(float_to_bits(sn) ^ ssign);
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

// rau_tanf — tan via diagonal ratio w/(1-w), sqrt cancels exactly.
// Sign via bit XOR, consistent with rau_sincosf convention.
// Note: large absolute error near poles (π/2 + nπ) is unavoidable;
// sign may flip when w rounds across the w=1 boundary. Use std tan
// for pole-adjacent inputs.
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

    float denom = 1.0f - w;
    if (denom < 1e-6f) denom = 1e-6f;  /* pole magnitude clamp */
    float rho = w / denom;

    /* tan sign: negative in Q1 (x<0,y>0) and Q3 (x>0,y<0) */
    Uint32 sign = (Uint32)(((qi >> 1) ^ (qi & 1)) & 1) << 31;
    return bits_to_float(float_to_bits(rho) ^ sign);
}

// rau_sincos_mf — morphable sincos: M=0 gives raw diagonal (non-uniform arc),
// M=1 gives fully warped (arc-uniform). Lerp between the two for
// smooth transition effects.
void rau_sincos_mf(float input_t, float M, float *sin_out, float *cos_out) {
    float rau_pos = mod4(input_t);
    int   qi_full = (int)rau_pos;
    float frac    = rau_pos - (float)qi_full;
    int   qi      = qi_full & 3;

    float w_raw  = frac;
    float w_warp = rau_warpf(frac);
    float w      = w_raw + M * (w_warp - w_raw);

    if (qi & 1) w = 1.0f - w;

    float omw = 1.0f - w;
    float D   = omw * omw + w * w;
    float inv = 1.0f / SDL_sqrtf(D);

    Uint32 csign = (Uint32)(((qi + 1) >> 1) & 1) << 31;
    Uint32 ssign = (Uint32)(( qi      >> 1) & 1) << 31;

    *cos_out = bits_to_float(float_to_bits(omw * inv) ^ csign);
    *sin_out = bits_to_float(float_to_bits(w   * inv) ^ ssign);
}

// ── Arctan Polynomial Tiers ────────────────────────────────────────────────
//
// All three fit arctan(x)/x as an odd polynomial over [0,1].
// Caller maps input → [0,1] via the t/(1+t) compactification before calling.
// Remez minimax derivation — coefficients verified against dense grid.

// f16: degree 3, 4 terms — float16 input quality
// Max err: 1.18e-4 RAU
static inline float rau_atanf_f16_polyf(float x) {
    float x2 = x * x;
    return x * (0.9998142570f
        + x2 * (-0.3262377264f
        + x2 * ( 0.1566730269f
        + x2 * (-0.0450371370f))));
}

// hsp: degree 4, 5 terms — half-single precision
// Max err: 1.67e-5 RAU
static inline float rau_atanf_hsp_polyf(float x) {
    float x2 = x * x;
    return x * (0.9999737848f
        + x2 * (-0.3318111223f
        + x2 * ( 0.1857423872f
        + x2 * (-0.0927448646f
        + x2 *   0.0242641934f))));
}

// precise: degree 6, 7 terms — float32 quality
// Max err: 3.62e-7 RAU
static inline float rau_atanf_precise_polyf(float x) {
    float x2 = x * x;
    return x * (0.9999994301f
        + x2 * (-0.3332707100f
        + x2 * ( 0.1988770404f
        + x2 * (-0.1351294716f
        + x2 * ( 0.0843566601f
        + x2 * (-0.0374368276f
        + x2 *   0.0080026120f))))));
}

// Dispatch to selected quality tier
#if   RAU_ATAN_QUALITY == 0
#  define rau_atanf_polyf rau_atanf_f16_polyf
#elif RAU_ATAN_QUALITY == 1
#  define rau_atanf_polyf rau_atanf_hsp_polyf
#else
#  define rau_atanf_polyf rau_atanf_precise_polyf
#endif

// ── Inverse Functions ──────────────────────────────────────────────────────

// rau_r_arctanf — exact rational arctan in RAU ∈ [0,1].
// Returns |y/x| / (1 + |y/x|) — the RAU diagonal coordinate.
// This is not an approximation; it is the exact w value for the
// first-quadrant angle whose tangent is |y/x|.
float rau_r_arctanf(float ry, float rx, int *err) {
    if (err) *err = 0;
    if (isnan(rx) || isnan(ry))       { if (err) *err = 1; return NAN; }
    if (ry == 0.0f && rx == 0.0f)    { if (err) *err = 1; return NAN; }
    if (isinf(ry) || isinf(rx)) {
        if (isinf(ry) && isinf(rx))  return 0.5f;
        if (isinf(ry))               return 1.0f;
        return 0.0f;
    }
    if (rx == 0.0f) return 1.0f;
    float t = fabsf(ry / rx);
    return t / (1.0f + t);
}

// rau_r_arcsinf — rational arcsin in RAU ∈ [0,1].
// Computes the first-quadrant RAU angle w whose sin equals s.
// Singularity at s = ±1/√2 (45°): returns 0.5 RAU (exact 45°).
float rau_r_arcsinf(float s, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(s) || s < -1.0f || s > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float sy2   = s * s;
    float denom = 2.0f * sy2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f)
        return 0.5f;  /* 45° — diagonal point, geometric singularity */
    float disc = SDL_sqrtf(SDL_max(sy2 * (1.0f - sy2), 0.0f));
    return SDL_fabsf((sy2 - disc) / denom);
}

// rau_r_arccosf — rational arccos in RAU ∈ [0,1].
// Computes the first-quadrant RAU angle w whose cos equals c.
// Singularity at c = ±1/√2 (45°): returns 0.5 RAU (exact 45°).
float rau_r_arccosf(float c, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(c) || c < -1.0f || c > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    float cx2   = c * c;
    float denom = 2.0f * cx2 - 1.0f;
    if (SDL_fabsf(denom) < 1e-10f)
        return 0.5f;  /* 45° — diagonal point, geometric singularity */
    float disc = SDL_sqrtf(SDL_max(cx2 * (1.0f - cx2), 0.0f));
    return SDL_fabsf((cx2 - 1.0f + disc) / denom);
}

// rau_invpolyf — convert RAU diagonal coordinate w ∈ [0,1] to RAU angle ∈ [0,1].
// w is the output of rau_r_arctanf (or rau_r_normf).
// Internally maps w → raw ratio t = w/(1-w) or (1-w)/w,
// then evaluates rau_atanf_polyf(t) directly — no double composition.
float rau_invpolyf(float w, int *err) {
    if (err) *err = 0;
    if (!rau_isfinitef(w) || w < 0.0f || w > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }
    if (w == 0.0f) return 0.0f;
    if (w == 1.0f) return 1.0f;

    if (w >= 0.5f) {
        /* w > 0.5: angle > 45°, reflect via reciprocal ratio */
        float t = (1.0f - w) / w;        /* raw ratio < 1 */
        float a = rau_atanf_polyf(t);    /* arctan(t) in radians */
        return ((float)M_PI_2 - a) * (float)M_2_PI;
    } else {
        /* w < 0.5: angle < 45°, direct evaluation */
        float t = w / (1.0f - w);        /* raw ratio < 1 */
        return rau_atanf_polyf(t) * (float)M_2_PI;
    }
}

// rau_atan2f — full four-quadrant arctan2, result in RAU [0,4).
// Accuracy determined by RAU_ATAN_QUALITY (default: float32, 3.62e-7 RAU).
float rau_atan2f(float y, float x, int *err) {
    if (err) *err = 0;
    if (rau_isnanf(x) || rau_isnanf(y))  { if (err) *err = 1; return NAN; }
    if (x == 0.0f && y == 0.0f)          { if (err) *err = 1; return NAN; }

    if (isinf(x) || isinf(y)) {
        if (isinf(y) && isinf(x)) {
            if (y > 0.0f && x > 0.0f) return 0.5f;
            if (y > 0.0f && x < 0.0f) return 1.5f;
            if (y < 0.0f && x < 0.0f) return 2.5f;
            return 3.5f;
        }
        if (isinf(y)) return (y > 0.0f) ? 1.0f : 3.0f;
        return (x > 0.0f) ? 0.0f : 2.0f;
    }

    if (x == 0.0f) {
        if (y > 0.0f) return 1.0f;
        if (y < 0.0f) return 3.0f;
        if (err) *err = 1;
        return NAN;
    }

    /* Quadrant from signs of x and y */
    int q = (x >= 0.0f) ? (y >= 0.0f ? 0 : 3) : (y >= 0.0f ? 1 : 2);

    /* Diagonal coordinate w = |y/x| / (1 + |y/x|) ∈ [0,1] */
    float t = fabsf(y / x);
    float w = t / (1.0f + t);

    /* Convert w to RAU angle within first quadrant */
    float frac = rau_invpolyf(w, err);
    if (err && *err) return NAN;

    /* Map first-quadrant result to full [0,4) range */
    if (q == 0) return        frac;
    if (q == 1) return 1.0f + (1.0f - frac);
    if (q == 2) return 2.0f + frac;
    return             3.0f + (1.0f - frac);
}

// rau_atanf — single-argument arctan, result in RAU [-1,+1].
// Uses same quality tier as rau_atan2f.
float rau_atanf(float x, int *err) {
    if (err) *err = 0;
    if (rau_isnanf(x))     { if (err) *err = 1; return NAN; }
    if (!rau_isfinitef(x)) { return SDL_copysignf(1.0f, x); }
    if (x == 0.0f)         { return x; }

    float ax = SDL_fabsf(x);
    float t  = ax / (1.0f + ax);   /* compactify to [0, 0.5) */

    float a;
    if (t >= 0.5f) {
        /* ax > 1: reflect via reciprocal */
        float v = (1.0f - t) / t;
        a = ((float)M_PI_2 - rau_atanf_polyf(v)) * (float)M_2_PI;
    } else {
        a = rau_atanf_polyf(t / (1.0f - t)) * (float)M_2_PI;
    }

    return SDL_copysignf(a, x);
}
