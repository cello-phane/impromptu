/* radicaltrig.h — Radical Angle Unit (RAU) trigonometry
 * Michael Ledesma (cello-phane@github)  |  Unlicense / Public Domain
 * Desmos: https://www.desmos.com/calculator/gkellct2v2
 *
 * 1 RAU = π/2 rad.  Full circle = 4 RAU.  Integer part = quadrant.
 * Diagonal ratio w = sin/(sin+cos) ∈ [0,1] parameterises one quadrant.
 *
 * degrees → RAU:  phi = deg  / 90.0f
 * radians → RAU:  phi = rad  * (2.0f / M_PI)
 *
 * RAU_ATAN_QUALITY (compile-time, default 2):
 *   0 = f16     ~1.18e-4 RAU   4 terms
 *   1 = hsp     ~1.67e-5 RAU   5 terms
 *   2 = precise ~3.62e-7 RAU   7 terms
 *
 * NON_UNIFORM_VEL (compile-time, default 0):
 *   0 = warp applied — arc-uniform, π/2 per RAU
 *   1 = raw diagonal  — no warp, non-uniform arc speed
 *   runtime blend: rau_sincos_mf(phi, M, s, c)  M∈[0,1]
 *
 * Accuracy:
 *   rau_sincosf       1.68e-7   float32
 *   rau_tanf          7.73e-6   float16  (pole sign unreliable near ±90°)
 *   rau_r_arctanf     ~1e-15    exact rational
 *   rau_r_arcsinf     ~1e-14    algebraic
 *   rau_r_arccosf     ~1e-14    algebraic
 *   rau_invpolyf      ~1.05e-4  approximate warp inverse
 *   rau_atan2f        ~3.62e-7  float32 (at default quality)
 *
 * Inverse functions accept int *err — set non-zero on invalid input, pass NULL to ignore.
 */

#ifndef RADICALTRIG_H
#define RADICALTRIG_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── §1 Forward ───────────────────────────────────────────────────────────── */

/* Warp polynomial: linear t ∈ [0,1] → arc-uniform w ∈ [0,1].
 * Remez minimax, leading coeff ≈ π/4. Max err: 1.61e-7. */
float rau_warpf(float t);

/* Simultaneous sin and cos from RAU angle. Single warp + single sqrt.
 * Prefer over separate rau_sinf/rau_cosf when both outputs are needed. */
void  rau_sincosf(float phi, float *sin_out, float *cos_out);

/* sin and cos from RAU angle — wrappers around rau_sincosf. */
float rau_sinf(float phi);
float rau_cosf(float phi);

/* tan from RAU angle. sqrt cancels in w/(1-w) — no sqrt required.
 * Pole clamp at ±90°/±270°: returns ±1e6, sign may be wrong near boundary. */
float rau_tanf(float phi);

/* Morphable sincos: M=0 raw diagonal, M=1 arc-uniform, M∈(0,1) blend. */
void  rau_sincos_mf(float phi, float M, float *sin_out, float *cos_out);

/* ── §2 Inverse ───────────────────────────────────────────────────────────── */

/* Exact rational arctan: w = |ry/rx| / (1 + |ry/rx|) ∈ [0,1].
 * Returns first-quadrant diagonal coordinate. No approximation. */
float rau_r_arctanf(float ry, float rx, int *err);

/* Algebraic arcsin from sin value → w ∈ [0,1].
 * Singularity at sin = ±1/√2 (45°): returns 0.5 exactly. */
float rau_r_arcsinf(float s, int *err);

/* Algebraic arccos from cos value → w ∈ [0,1].
 * Singularity at cos = ±1/√2 (45°): returns 0.5 exactly. */
float rau_r_arccosf(float c, int *err);

/* Approximate inverse of rau_warpf: w ∈ [0,1] → t ∈ [0,1).
 * Accuracy set by RAU_ATAN_QUALITY. Max err at default: ~3.62e-7 RAU. */
float rau_invpolyf(float w, int *err);

/* Four-quadrant phase recovery → RAU ∈ [0,4). Analogous to atan2(y,x).
 * 0=+x, 1=+y, 2=-x, 3=-y. NaN input → err=1. */
float rau_atan2f(float y, float x, int *err);

/* Principal-value arctan → RAU ∈ [-1,+1]. Analogous to atan(x)/(π/2).
 * ±inf → ±1. NaN → err=1. Preserves signed zero. */
float rau_atanf(float x, int *err);

/* ── §3 Conversion helpers ────────────────────────────────────────────────── */

/* Wrap unsigned atan2 result [0,4) to signed (−π, +π] radians. */
float rau_atan2_signed_radians(float phi_rau);

/* Wrap unsigned atan2 result [0,4) to degrees [0°, 360°). */
float rau_atan2_signed_degs(float phi_rau);

#ifdef __cplusplus
}
#endif

#endif /* RADICALTRIG_H */

/* LICENSE — Unlicense / Public Domain
 * This is free and unencumbered software released into the public domain.
 * See <https://unlicense.org> for full terms.
 * Author: Michael Ledesma (cello-phane@github)
 */
