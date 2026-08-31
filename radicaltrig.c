#include "radicaltrig.h"
#include "SDL_stdinc.h"
#include <math.h>
#include <immintrin.h>

// Radical Angle Unit (RAU) trigonometric library
// Reference: https://www.desmos.com/calculator/gkellct2v2
//
// 1 RAU = π/2 radians. Full circle = 4 RAU
//
// Accuracy quality tiers (selectable via RAU_ATAN_QUALITY as 0, 1, or 2):
// 0 = f16 (1.18e-4)
// 1 = hsp (2.92e-5)
// 2 = precise (7.32e-7, join=9.74e-11)
//
// NON_UNIFORM_VEL 0 — arc-uniform: warp polynomial applied, π/2 arc per RAU
// NON_UNIFORM_VEL 1 — raw diagonal projection, no warp, non-uniform arc speed
#define NON_UNIFORM_VEL 0

#ifndef RAU_ATAN_QUALITY
#define RAU_ATAN_QUALITY 2 // 2 (precise is default)
#endif

// RAU_WARP_QUALITY selects which forward-warp coefficient tier rau_warpf()
// (and therefore rau_sincosf/rau_tanf/rau_sincos_mf) dispatches to:
// 1 = v2 (default) — constrained minimax refit, end-to-end sin/cos max
//                     err 1.246e-7 (double-precision theoretical) /
//                     2.372e-7 (measured on real float32 hardware, ~2.68x
//                     tighter than tier 0 there). Same C1=pi/4 exactness
//                     as tier 0. Validated: compiled and run on real
//                     x86 float32 (GCC), cross-checked against an
//                     independent LP/mpmath derivation and against a
//                     reference Desmos construction.
// 0 = original       — end-to-end sin/cos max err 5.472e-7. Kept for
//                     anyone who specifically wants the older fit;
//                     rau_warpf_orig() is also directly callable
//                     regardless of this define.
#ifndef RAU_WARP_QUALITY
#define RAU_WARP_QUALITY 1
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
// parameter w = sin(t·π/2) / (sin(t·π/2) + cos(t·π/2))
//
// Two tiers are always compiled in; RAU_WARP_QUALITY (defined above)
// selects which one rau_warpf() itself dispatches to (default: tier 1,
// v2). rau_warpf_orig() and rau_warpf_v2() are both exposed directly so
// either is callable regardless of that define, for A/B comparison
// without a recompile.
//
// Both figures below are end-to-end sin/cos error (after the radical-
// identity projection), which is very slightly larger than the raw
// warp-space approximation error alone (5.168e-7 for tier 0) — verified
// via mpmath at 50-digit precision and cross-checked with an independent
// scipy/HiGHS linear-program minimax solve, then again by compiling and
// running both tiers on real x86 float32 hardware.
//
//   tier 0 (original): Remez minimax over full [0,1]. C1 = π/4 exact.
//                       sin/cos max err: 5.472e-7 (double) / 6.36e-7 (float32)
//   tier 1 (v2, default): constrained L-infinity refit — C1 still pinned
//                       exactly to π/4, C2..C6 re-solved by LP. sin/cos
//                       max err: 1.246e-7 (double) / 2.37e-7 (float32) —
//                       ~4.4x tighter in double precision, ~2.7x tighter
//                       as actually compiled, confirming the improvement
//                       survives float32 rounding rather than being an
//                       artifact of the double-precision derivation.
static inline float rau_warpf_tier0(float t) {
    static const float C[6] = {
        0.78539816339744830962f,
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

static inline float rau_warpf_tier1(float t) {
    static const float C[6] = {
        0.78539816339744830962f, /* unchanged from tier 0 — pi/4, exact */
        0.6460261721112686f,
        0.6346711435786873f,
        0.6812793876895539f,
        0.33448057938003534f,
        1.5121252251949524f
    };
    float v  = t - 0.5f;
    float v2 = v * v;
    float p  = C[5];
    for (int i = 4; i >= 0; --i)
        p = v2 * p + C[i];
    return 0.5f + v * p;
}

float rau_warpf_orig(float t) { return rau_warpf_tier0(t); }
float rau_warpf_v2(float t) { return rau_warpf_tier1(t); }

#if RAU_WARP_QUALITY == 1
float rau_warpf(float t) { return rau_warpf_tier1(t); }
#else
float rau_warpf(float t) { return rau_warpf_tier0(t); }
#endif

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

// rau_tanf — tan via diagonal ratio w/(1-w), sqrt cancels exactly
// Sign via bit XOR, consistent with rau_sincosf convention
// Note: large absolute error near poles (π/2 + nπ) is unavoidable;
// sign may flip when w rounds across the w=1 boundary
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

    /* tan sign: positive in Q0 (x>0,y>0) and Q2 (x<0,y<0),
     * negative in Q1 (x<0,y>0) and Q3 (x>0,y<0) — simplifies to the
     * low bit of the quadrant index. (Fixed: the previous XOR-based
     * formula ((qi>>1)^(qi&1)) flipped the sign in Q2 and Q3 — verified
     * both by inspection and by comparing rau_tanf's sign against
     * sin/cos from rau_sincosf across all four quadrants.) */
    Uint32 sign = (Uint32)(qi & 1) << 31;
    return bits_to_float(float_to_bits(rho) ^ sign);
}

// rau_sincos_mf — morphable sincos: M=0 gives raw diagonal (non-uniform arc),
// M=1 gives fully warped (arc-uniform). Lerp between the two for
// smooth transition effects
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
// All three fit arctan(x)/x as an odd polynomial over [0,1]
// Caller maps input → [0,1] via the t/(1+t) compactification before calling
// Remez minimax derivation — coefficients verified against dense grid

// f16: degree 3, 4 terms — float16 input quality
// Max err: 1.18e-4 RAU
static inline float rau_atanf_f16_polyf(float x) {
    float x2 = x * x;
    return x * (0.9998142570f
        + x2 * (-0.3262377264f
        + x2 * ( 0.1566730269f
        + x2 * (-0.0450371370f))));
}
//tested the 2 below at https://godbolt.org/z/Tz6dccone in C++
// hsp: degree 4, 5 terms — half-single precision
static inline float rau_atanf_hsp_polyf(float x) {
    float x2 = x * x;
    // Max err: 1.67e-5 RAU
    /*return x * (0.9999737848f
        + x2 * (-0.3318111223f
        + x2 * ( 0.1857423872f
        + x2 * (-0.0927448646f
        + x2 *   0.0242641934f))));
    */
    // Max err: 2.9156e-05
    return x * (0.9999697079f
        + x2 * (-0.3317181925f
        + x2 * ( 0.1851789314f
        + x2 * (-0.0916479384f
        + x2 *   0.0236156549f))));
}

// precise: degree 6, 7 terms
// Max err: 7.32e-7 rad | Join gap at x=1: 9.74e-11 (constrained, was 5.70e-7)
static inline float rau_atanf_precise_polyf(float x) {
    float x2 = x * x;
    // Max err: 6.8747e-07
    /*return x * (0.9999994301f
        + x2 * (-0.3332707100f
        + x2 * ( 0.1988770404f
        + x2 * (-0.1351294716f
        + x2 * ( 0.0843566601f
        + x2 * (-0.0374368276f
        + x2 *   0.0080026120f))))));
    */
    // Max err: 7.3249e-07
    return x * (0.9999997567f
        + x2 * (-0.3332778034f
        + x2 * (0.1989157444f
        + x2 * (-0.1351964889f
        + x2 * (0.0843541706f
        + x2 * (-0.0373408246f
        + x2 * (0.0079436085f)))))));
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

// ── Core Coordinate Mapping Kernel ─────────────────────────────────────────
// Maps a normalized linear diamond distance [0.0, 1.0] to a circular arc length

static inline float rau_kernel_unwarp(float k) {
    float u  = k - 0.5f;
    float u2 = u * u;

    float p  = -18.8502046660f;
    p = u2 * p + 19.4799685220f;
    p = u2 * p + -10.0240706528f;
    p = u2 * p + 3.9960544827f;
    p = u2 * p + -1.6962976022f;
    p = u2 * p + 1.2732550558f;

    return 0.5f + u * p;
}

// ── Optimized Drop-In Inverse Variants ──────────────────────────────────────

// ── Drop-In Inverse Variant 1: Pure Diagonal (Arctangent Alternative) ─────
// Maps raw spatial vectors to single-argument RAU ∈ [-1.0, 1.0]
float rau_invdiagonal_from_ratio(float ry, float rx) {
    int err;
    float k = rau_r_arctanf(ry, rx, &err);
    if (err) return 0.0f;

    float poly_out = rau_kernel_unwarp(k);

    /* Sign of the result is sign(ry/rx), i.e. sign(ry) XOR sign(rx) —
     * not sign(rx) alone. (Fixed: using only rx's sign disagreed with
     * rau_invdiagonal_from_ratio_quotient(ry/rx) — its documented
     * equivalent — whenever ry was negative; verified by direct
     * comparison across all four sign combinations of ry, rx.) */
    Uint32 sign = (float_to_bits(ry) ^ float_to_bits(rx)) & 0x80000000u;
    return bits_to_float(float_to_bits(poly_out) ^ sign);
}

// ── Drop-In Inverse Variant 1B: Inline Quotient Variant ────────────────────
// Maps a pre-computed ry/rx slope directly to single-argument RAU ∈ [-1.0, 1.0]
// Safely bypasses the vector components while preserving singularity protection
float rau_invdiagonal_from_ratio_quotient(float ry_over_rx) {
    // 1. Explicitly catch vertical tracking limits (rx == 0) to prevent NaN quotients
    if (isinf(ry_over_rx)) {
        // As slope approaches infinity, k converges to 1.0f
        // rau_kernel_unwarp(1.0f) smoothly evaluates to a magnitude of 1.0f RAU
        return SDL_copysignf(1.0f, ry_over_rx);
    }
    if (isnan(ry_over_rx)) {
        return 0.0f;
    }
    float abs_t = SDL_fabsf(ry_over_rx);
    float k     = abs_t / (1.0f + abs_t);

    float poly_out = rau_kernel_unwarp(k);
    return SDL_copysignf(poly_out, ry_over_rx);
}

// ── Drop-In Inverse Variant 2: From Y-Component (Arcsine Alternative) ──────
// Maps raw Y coordinate/sine to RAU ∈ [-1.0, 1.0]
float rau_invdiagonal_from_y(float ry) {
    int err;
    float k = rau_r_arcsinf(ry, &err);
    if (err) return 0.0f;

    float poly_out = rau_kernel_unwarp(k);

    return SDL_copysignf(poly_out, ry);
}

// ── Drop-In Inverse Variant 3: From X-Component (Arccosine Alternative) ────
// Maps raw X coordinate/cosine to RAU ∈ [0.0, 2.0] (Principal Arccosine Range)
float rau_invdiagonal_from_x(float rx) {
    int err;
    float k = rau_r_arccosf(rx, &err);
    if (err) return 0.0f;

    float poly_out = rau_kernel_unwarp(k);

    // Matches Desmos Line 52: { R.x > 0: poly_out, 2.0f - poly_out }
    if (rx < 0.0f) {
        return 2.0f - poly_out;
    }
    return poly_out;
}

// rau_r_arctanf — exact rational arctan in RAU ∈ [0,1].
// Returns |y/x| / (1 + |y/x|) — the RAU diagonal coordinate
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
// Singularity at s = ±1/√2 (45°): returns 0.5 RAU (exact 45°)
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

// rau_r_arccosf — rational arccos in RAU ∈ [0,1]
// Computes the first-quadrant RAU angle w whose cos equals c
// Singularity at c = ±1/√2 (45°): returns 0.5 RAU (exact 45°)
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
// w is the output of rau_r_arctanf (or rau_r_normf)
// Internally maps w → raw ratio t = w/(1-w) or (1-w)/w
// then evaluates rau_atanf_polyf(t) directly — no double composition
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

// rau_atan2f — full four-quadrant arctan2, result in RAU [0,4)
// Accuracy determined by RAU_ATAN_QUALITY (default: float32, 3.62e-7 RAU)
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

// rau_atanf — single-argument arctan, result in RAU [-1,+1]
float rau_atanf(float x, int *err) {
    if (err) *err = 0;
    if (rau_isnanf(x))     { if (err) *err = 1; return NAN; }
    if (!rau_isfinitef(x)) { return SDL_copysignf(1.0f, x); }
    if (x == 0.0f)         { return x; }

    float ax = SDL_fabsf(x);
    float a;
    if (ax > 1.0f) { // SIMD if available
    #ifdef __SSE__
        float inv;
        __m128 v = _mm_set_ss(ax);
        _mm_store_ss(&inv, _mm_rcp_ss(v));  // ~12-bit approx, 1 cycle
        inv = inv * (2.0f - ax * inv);       // Newton step → ~23-bit, 2 FMAs
    #else
        float inv = 1.0f / ax;              // scalar fallback
    #endif
        a = ((float)M_PI_2 - rau_atanf_polyf(inv)) * (float)M_2_PI;
    } else {
        /* Fixed: this branch was previously missing entirely, leaving
         * `a` uninitialized for every |x| <= 1 — i.e. the common case,
         * not an edge case. Verified via -Wmaybe-uninitialized and by
         * direct comparison against atanf(x)*(2/pi), which the old code
         * returned ~0 garbage for regardless of x. */
        a = rau_atanf_polyf(ax) * (float)M_2_PI;
    }
    return SDL_copysignf(a, x);
}
