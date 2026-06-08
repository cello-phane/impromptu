#include "radicaltrig.h"
#include "SDL_stdinc.h"
#include <math.h>

// Graph of sin/cos/tan and inverse functions -
// https://www.desmos.com/calculator/gkellct2v2

// NON_UNIFORM_VEL 0 — arc-uniform: warp polynomial applied, constant arc speed
// (π/2 per unit parameter) NON_UNIFORM_VEL 1 — raw linear diagonal projection,
// no warp, fast, non-uniform in arc
#define NON_UNIFORM_VEL 0

// ── Helper Functions ────────────────────────────────────────
static float mod4(float a) {
  float r = SDL_fmodf(a, 4.0f);
  if (r < 0.0f)
    r += 4.0f;
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

static inline int rau_isfinitef(float x) { return isfinite(x); }
static inline int rau_isnanf(float x) { return isnan(x); }
static inline int rau_in_unit_range(float x) { return x >= -1.0f && x <= 1.0f; }

static inline float rau_clampf(float x, float lo, float hi) {
  return x < lo ? lo : (x > hi ? hi : x);
}

static inline float rau_sanitize_unit(float x) {
  if (x < -1.0f)
    return -1.0f;
  if (x > 1.0f)
    return 1.0f;
  return x;
}

// ── Convert range of unsigned atan2 result to signed (-pi,+pi] ─────────────
float rau_atan2_signed_radians(float phi_rau) {
  float wrapped = fmodf(phi_rau + 2.0f, 4.0f);
  if (wrapped < 0.0f)
    wrapped += 4.0f;
  wrapped -= 2.0f;
  return wrapped * (float)M_PI_2;
}

float rau_atan2_signed_degs(float phi_rau) {
  float deg = fmodf(phi_rau * 90.0f, 360.0f);
  if (deg < 0.0f)
    deg += 360.0f;
  return deg;
}

// ── Forward Trigonometric Functions ────────────────────────────────────────
// float rau_warpf(float t) {
//   static const float C[6] = {
//       0.78539732f, /* ≈ π/4 — Remez-confirmed leading coefficient */
//       0.64607089f, 0.63401085f, 0.68518412f, 0.32482231f, 1.52006419f};
//   float v = t - 0.5f; /* centre on 0 to exploit odd symmetry */
//   float v2 = v * v;
//   float p = C[5]; /* Nested sum of products in v² (even powers only) */
//   for (int i = 4; i >= 0; --i)
//     p = v2 * p + C[i];
//   return 0.5f + v * p;
// }

// Wider bit version than above(but slightly less precise):
float rau_warpf(float t)
{
    static const float C[6] = {
        0.78539816339744830962f, /* ≈ π/4 — Remez-confirmed leading coefficient */
        0.64607158024987317298f, 0.63401589172451679138f,
        0.68515350354689586789f, 0.32501622369042378935f,
        1.51901679307446258196f};
    float v  = t - 0.5f; /* centre on 0 to exploit odd symmetry */
    float v2 = v * v;
    float p = C[5]; /* Nested sum of products in v² (even powers only) */
    for (int i = 4; i >= 0; --i)
        p = v2 * p + C[i];
    return 0.5f + v * p;
}

void rau_sincosf(float input_t, float *sin_out, float *cos_out) {
  /* Range reduction: fold into [0,4), extract integer quadrant */
  float rau_pos = mod4(input_t);
  int qi_full = (int)rau_pos;
  float frac = rau_pos - (float)qi_full; /* fractional part ∈ [0,1) */
  int qi = qi_full & 3;                  /* quadrant: 0,1,2,3 */

/* Arc-length correction (compile-time switch) */
#if NON_UNIFORM_VEL
  float w = frac; /* raw diagonal: non-uniform arc speed */
#else
  float w = rau_warpf(frac); /* warped: uniform arc speed = π/2 per RAU */
#endif

  /* Odd-quadrant reversal: w increases [0→1] in all quadrants after this */
  if (qi & 1)
    w = 1.0f - w;

  /* Diagonal point to unit circle — single sqrt, both outputs */
  float omw = 1.0f - w;
  float D =
      omw * omw + w * w; /* squared distance from origin to diagonal point */
  float inv = 1.0f / SDL_sqrtf(D); /* normalisation factor */
  float cs = omw * inv;            /* cos = (1-w)/sqrt(D) */
  float sn = w * inv;              /* sin = w /sqrt(D) */

  Uint32 csign = (Uint32)(((qi + 1) >> 1) & 1) << 31;
  Uint32 ssign = (Uint32)((qi >> 1) & 1) << 31;

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
  int qi_full = (int)rau_pos;
  float frac = rau_pos - (float)qi_full;
  int qi = qi_full & 3;

#if NON_UNIFORM_VEL
  float w = frac;
#else
  float w = rau_warpf(frac);
#endif

  if (qi & 1)
    w = 1.0f - w;

  float denom = 1.0f - w; /* = cos (without normalisation) */
  if (denom < 1e-6f)
    denom = 1e-6f; /* pole clamp — see note above */

  float t = w / denom; /* = tan (sqrt cancels) */
  if (qi & 1)
    t = -t; /* sign: tan negative in Q1,Q3 */
  return t;
}

void rau_sincos_mf(float input_t, float M, float *sin_out, float *cos_out) {
  float rau_pos = mod4(input_t);
  int qi_full = (int)rau_pos;
  float frac = rau_pos - (float)qi_full;
  int qi = qi_full & 3;

  /* Lerp between raw (M=0) and warped (M=1) */
  float w_raw = frac;
  float w_warp = rau_warpf(frac);
  float w = w_raw + M * (w_warp - w_raw);

  if (qi & 1)
    w = 1.0f - w;

  float omw = 1.0f - w;
  float D = omw * omw + w * w;
  float inv = 1.0f / SDL_sqrtf(D);
  float cs = omw * inv;
  float sn = w * inv;

  Uint32 csign = (Uint32)(((qi + 1) >> 1) & 1) << 31;
  Uint32 ssign = (Uint32)((qi >> 1) & 1) << 31;

  *cos_out = bits_to_float(float_to_bits(cs) ^ csign);
  *sin_out = bits_to_float(float_to_bits(sn) ^ ssign);
}

float rau_r_arctanf(float ry, float rx, int *err) {
  if (err)
    *err = 0;
  if (isnan(rx) || isnan(ry)) {
    if (err)
      *err = 1;
    return NAN;
  }
  if (ry == 0.0f && rx == 0.0f) {
    if (err)
      *err = 1;
    return NAN;
  }
  if (isinf(ry) || isinf(rx)) {
    if (isinf(ry) && isinf(rx))
      return 0.5f;
    if (isinf(ry))
      return 1.0f;
    return 0.0f;
  }
  if (rx == 0.0f)
    return 1.0f;

  float t = fabsf(ry / rx);
  return t / (1.0f + t);
}

float rau_r_arcsinf(float s, int *err) {
  if (err)
    *err = 0;
  if (!rau_isfinitef(s) || s < -1.0f || s > 1.0f) {
    if (err)
      *err = 1;
    return NAN;
  }
  float sy2 = s * s;
  float denom = 2.0f * sy2 - 1.0f;
  if (SDL_fabsf(denom) < 1e-10f)
    return 0.5f; /* 45° singularity guard */
  float disc = SDL_sqrtf(SDL_max(sy2 * (1.0f - sy2), 0.0f));
  return SDL_fabsf((sy2 - disc) / denom);
}

float rau_r_arccosf(float c, int *err) {
  if (err)
    *err = 0;
  if (!rau_isfinitef(c) || c < -1.0f || c > 1.0f) {
    if (err)
      *err = 1;
    return NAN;
  }
  float cx2 = c * c;
  float denom = 2.0f * cx2 - 1.0f;
  if (SDL_fabsf(denom) < 1e-10f)
    return 0.5f; /* 45° singularity guard */
  float disc = SDL_sqrtf(SDL_max(cx2 * (1.0f - cx2), 0.0f));
  return SDL_fabsf((cx2 - 1.0f + disc) / denom);
}

// Remez minimax, arctan(x)/x over [0,1], degree 3 — 4 terms
// Max err: 1.18e-4 RAU — float16 quality, minimal cost
static inline float rau_atanf_sp_polyf(float x)
{
    float x2 = x * x;
    return x * (0.9998142570f
        + x2 * (-0.3262377264f
        + x2 * ( 0.1566730269f
        + x2 * (-0.0450371370f))));
}

// Remez minimax, arctan(x)/x over [0,1], degree 4 — 5 terms
// Max err: 1.67e-5 RAU — half-single precision quality
static inline float rau_atanf_hsp_polyf(float x)
{
    float x2 = x * x;
    return x * (0.9999737848f
        + x2 * (-0.3318111223f
        + x2 * ( 0.1857423872f
        + x2 * (-0.0927448646f
        + x2 *   0.0242641934f))));
}

static inline float rau_atanf_polyf(float x)
{
    float x2 = x * x;
    return x * (0.9999994300813833f
        + x2 * (-0.3332707100403806f
        + x2 * ( 0.1988770404319326f
        + x2 * (-0.1351294715695766f
        + x2 * ( 0.0843566600833229f
        + x2 * (-0.0374368276347488f
        + x2 *   0.0080026119641338f))))));
}

static inline int rau_quadrantf(float y, float x)
{
    if (x >= 0.0f) return (y >= 0.0f) ? 0 : 3;
    return (y >= 0.0f) ? 1 : 2;
}

static inline float rau_r_normf(float y, float x)
{
    if (x == 0.0f) return 1.0f;
    float t = fabsf(y / x);
    return t / (1.0f + t);
}

/* Reduced-domain arctan helper, safe only for finite x/y when x != 0. */
float rau_invpolyf(float w, int *err)
{
    if (err) *err = 0;

    if (!rau_isfinitef(w) || w < 0.0f || w > 1.0f) {
        if (err) *err = 1;
        return NAN;
    }

    if (w == 0.0f) return 0.0f;
    if (w == 1.0f) return 1.0f;

    if (w >= 0.5f) {
        float u = (1.0f - w) / w;
        float a = rau_atanf_polyf(u);
        return (M_PI_2 - a) * M_2_PI;
    } else {
        float u = w / (1.0f - w);
        return rau_atanf_polyf(u) * (2.0f / M_PI);
    }
}

/* Returns RAU angle in [0,4). 1 RAU = pi/2 radians. */
float rau_atan2f(float y, float x, int *err)
{
    if (err) *err = 0;
    if (rau_isnanf(x) || rau_isnanf(y)) {
        if (err) *err = 1;
        return NAN;
    }
    if (x == 0.0f && y == 0.0f) {
        if (err) *err = 1;
        return NAN;
    }

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

    int q = rau_quadrantf(y, x);
    float w = rau_r_normf(y, x);
    float t = rau_invpolyf(w, err);
    if (err && *err) return NAN;

    if (q == 0) return 0.0f + t;
    if (q == 1) return 1.0f + (1.0f - t);
    if (q == 2) return 2.0f + t;
    return 3.0f + (1.0f - t);
}

/* Returns a signed normalized RAU angle in [-1, 1]. */
float rau_atanf(float x, int *err)
{
    if (err) *err = 0;
    if (rau_isnanf(x)) {
        if (err) *err = 1;
        return NAN;
    }
    if (!rau_isfinitef(x)) {
        return copysignf(1.0f, x);
    }
    if (x == 0.0f) {
        return x;
    }

    float ax = SDL_fabsf(x);
    float u = ax / (1.0f + ax);

    float a;
    if (u >= 0.5f) {
        float v = (1.0f - u) / u;
        a = (M_PI_2 - rau_atanf_polyf(v)) * (2.0f / M_PI);
    } else {
        float v = u / (1.0f - u);
        a = rau_atanf_polyf(v) * (2.0f / M_PI);
    }

    return SDL_copysignf(a, x);
}
