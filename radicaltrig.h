/* radicaltrig.h
 * =============
 * Radical Trigonometry
 * Michael Ledesma (cello-phane@github)  |  Unlicense / Public Domain
 *
 * Desmos reference: https://www.desmos.com/calculator/gkellct2v2
 *
 * A trigonometry library built on the diagonal ratio
 * parameterisation of the unit circle. All functions use only polynomial
 * arithmetic, one sqrt (which cancels in tan), and IEEE 754 bit operations.
 *
 * ── Radical Angle Unit (RAU) ────────────────────────────────────────────────
 *
 *   0 RAU = 0°     1 RAU = 90°     2 RAU = 180°     4 RAU = 360°
 *
 *   The diagonal ratio w = sin/(sin+cos) ∈ [0,1] parameterises one quadrant.
 *   Integer part of the RAU angle is the quadrant index (0–3).
 *   Fractional part is the within-quadrant parameter t ∈ [0,1).
 *
 * ── Input conversion ────────────────────────────────────────────────────────
 *
 *   degrees → RAU:  phi = degrees * (1.0f / 90.0f)
 *   radians → RAU:  phi = radians * (2.0f / M_PI)
 *
 *   The radian conversion is a presentation conversion; forward trig still
 *   operates on RAU phase. When NON_UNIFORM_VEL=1, the internal phase-to-arc
 *   mapping uses the raw linear diagonal projection instead of the warp.
 *
 * ── NON_UNIFORM_VEL compile toggle ──────────────────────────────────────────
 *
 *   #define NON_UNIFORM_VEL 0  (default)
 *     Warp polynomial applied. Arc-uniform velocity.
 *     Correct for FFT twiddles, motor control, camera pans, PLL phase detection.
 *
 *   #define NON_UNIFORM_VEL 1
 *     Raw linear diagonal projection. No warp.
 *     Faster, but arc speed is non-uniform.
 *
 *   For runtime blend between modes, use rau_sincos_mf(phi, M, s, c)
 *   where M=0 gives NON_UNIFORM_VEL=1 behaviour and M=1 gives NON_UNIFORM_VEL=0.
 *
 * ── Accuracy ────────────────────────────────────────────────────────────────
 *
 *   Function          Max error   Notes
 *   rau_sincosf       ~1.68e-7    float32 quality
 *   rau_tanf          ~7.73e-6    limited by warp polynomial
 *   rau_r_arctanf     ~6.73e-16    machine-epsilon-level rounding only
 *   rau_r_arcsinf     ~1.13e-14    algebraic inverse + rounding
 *   rau_r_arccosf     ~1.13e-14    algebraic inverse + rounding
 *   rau_invpolyf      ~1.05e-4    approximate inverse of rau_warpf
 *   rau_atan2f        ~1.01e-4    full RAU phase recovery
 *
 * ── Error codes ─────────────────────────────────────────────────────────────
 *
 *   Inverse functions accept an optional int *err pointer.
 *   Pass NULL to ignore. Non-zero on invalid input (NaN, out-of-range, etc.)
 *
 * ── References ──────────────────────────────────────────────────────────────
 *
 *   Origin of diagonal parameterisation: John Gabriel
 *   Remez algorithm: Chebyshev equioscillation theorem
 *   IEEE 754 sign-bit XOR: avoids branches on hot path
 *   atan2 is four-quadrant inverse tangent; RAU phase recovery is unsigned
 */

#ifndef RADICALTRIG_H
#define RADICALTRIG_H

#ifdef __cplusplus
extern "C" {
#endif

/* ── Compile-time configuration ─────────────────────────────────────────── */

/* NON_UNIFORM_VEL: see file header for full explanation.
 * Define before including this header to override, or set via compiler flag:
 *   -DNON_UNIFORM_VEL=1
 */
#ifndef NON_UNIFORM_VEL
#  define NON_UNIFORM_VEL 0
#endif

/* ══════════════════════════════════════════════════════════════════════════
 * §1  FORWARD FUNCTIONS
 * ══════════════════════════════════════════════════════════════════════════ */

/* rau_warpf — warp polynomial, maps linear t ∈ [0,1] to arc-length w ∈ [0,1].
 *
 * Approximates the exact diagonal-to-arc conversion:
 *   w_exact(t) = sin(t·π/2) / (sin(t·π/2) + cos(t·π/2))
 *
 * Without this correction the diagonal parameterisation produces non-uniform
 * arc speed — 2× faster at 45° than at 0° or 90°.  Applying warp makes arc
 * speed exactly π/2 per RAU throughout the quadrant.
 *
 * Implementation: Nested sum-of-products evaluation in v² where v = t − 0.5.
 *   v² rather than t exploits the odd symmetry of (w − 0.5) around t = 0.5,
 *   halving the effective polynomial degree for the same accuracy.
 *
 * Coefficients: Remez minimax.  WC[0] ≈ π/4 — independently confirmed by the
 *   Remez algorithm, matching the analytically correct leading coefficient of
 *   the Taylor series for w_exact(t) around t = 0.5.
 *
 * Fixed point:    warp(0.5) = 0.5  exactly (by odd symmetry of the polynomial)
 * Boundary:       warp(0) ≈ 0,  warp(1) ≈ 1
 * Monotone:       yes, guaranteed by coefficient signs
 * Max error:      1.61e-7  (~22.6 effective bits, float32 quality)
 */
float rau_warpf(float t);

/* rau_sincosf — simultaneous cos and sin from RAU angle.
 *
 * Computes both outputs in a single call — cheaper than two separate calls
 * because both share the same warp evaluation and the same single sqrt.
 * No libm transcendental functions required.
 *
 * Algorithm:
 *   1. Range reduction: mod4 folds any input to [0,4), extracts quadrant q ∈ {0,1,2,3}
 *   2. Warp (NON_UNIFORM_VEL=0): maps fractional t → arc-length w via rau_warpf
 *   3. Odd-quadrant reversal: Q1,Q3 flip w → 1−w so w increases [0→1] in all quadrants
 *   4. Diagonal to circle — single sqrt for both outputs:
 *        D   = (1−w)² + w²  =  1 − 2w(1−w),   D ∈ [0.5, 1]
 *        cos = (1−w) / sqrt(D)
 *        sin =    w  / sqrt(D)
 *   5. Sign correction via IEEE 754 sign-bit XOR (zero-cost, branch-free):
 *        cos negative in Q1,Q2:  (qi+1)>>1 & 1  →  sign pattern [0,1,1,0]
 *        sin negative in Q2,Q3:   qi>>1    & 1  →  sign pattern [0,0,1,1]
 *
 * The type-punning for XOR uses SDL_memcpy / memcpy, which is defined behaviour
 * in C99 and C11.  In C++, prefer std::bit_cast (C++20) or memcpy.
 *
 * Input:         RAU — 0 to 4 = full revolution
 * *cos_out:      cos of the angle
 * *sin_out:      sin of the angle
 * Pythagorean:   cos²+sin² = 1 to within 1e-6 (exact at cardinal points)
 * Max error:     1.68e-7  (~22.5 effective bits, float32 quality)
 *
 * If both outputs are needed, always prefer rau_sincosf over calling
 * rau_sinf + rau_cosf separately — that would compute warp and sqrt twice.
 */
void rau_sincosf(float input_t, float *sin_out, float *cos_out);

/* rau_sinf — sin from RAU angle.
 *
 * Wrapper around rau_sincosf.  If both sin and cos are needed at the same
 * angle, call rau_sincosf directly — it produces both from a single sqrt.
 *
 * Input:     RAU (0 to 4 = full revolution)
 * Max error: 1.68e-7  (float32 quality)
 */
float rau_sinf(float x);

/* rau_cosf — cos from RAU angle.
 *
 * Wrapper around rau_sincosf.  If both sin and cos are needed at the same
 * angle, call rau_sincosf directly — it produces both from a single sqrt.
 *
 * Input:     RAU (0 to 4 = full revolution)
 * Max error: 1.68e-7  (float32 quality)
 */
float rau_cosf(float x);

/* rau_tanf — tan from RAU angle.  No sqrt required.
 *
 * tan = sin/cos = (w/sqrt(D)) / ((1−w)/sqrt(D)) = w/(1−w)
 *
 * The sqrt(D) in numerator and denominator cancels exactly, making rau_tanf
 * strictly cheaper than rau_sincosf: one fewer sqrt call.
 *
 * The ratio ρ = w/(1−w) is the odds ratio, which is also the phase variable
 * in the transcendental-free PLL phase detector:
 *   ρ_e = (ρ_r − ρ_v) / (1 + ρ_r·ρ_v)  =  tan(φ_r − φ_v)
 * This is a Möbius transformation — exact rational arithmetic on the ρ values.
 *
 * Pole behaviour at 90° / 270°:
 *   The clamp `denom < 1e-6 → 1e-6` returns ±1e6 instead of ±inf.
 *   For signal processing where the sign of approach to the pole matters,
 *   replace the clamp in the .c source with:
 *     denom = SDL_copysignf(SDL_max(SDL_fabsf(denom), 1e-6f), denom);
 *
 * Note: radian conversion requires NON_UNIFORM_VEL=0 for correctness.
 *   NON_UNIFORM_VEL=1 returns tan of the raw diagonal angle, not the arc angle.
 *   These differ by the warp correction and are not interchangeable.
 *
 * Input:     RAU (0 to 4 = full revolution)
 * Max error: ~7.73e-6  (float16 quality — limited by warp polynomial)
 */
float rau_tanf(float x);

/* rau_sincos_mf — sin/cos with runtime warp blend parameter M.
 *
 * Generalises the compile-time NON_UNIFORM_VEL toggle to a runtime float:
 *   M = 0.0f  →  raw linear diagonal  (NON_UNIFORM_VEL=1 behaviour)
 *   M = 1.0f  →  warp polynomial      (NON_UNIFORM_VEL=0 behaviour)
 *   M ∈ (0,1) →  smooth lerp between the two — continuously variable transition
 *
 * Typical uses:
 *   Soft motor landing:   M=1 during transit (precise arc speed), ramp to M=0
 *                         near the target (natural deceleration cushion)
 *   Animation easing:     M controls tightness of the arc trajectory
 *   CVT blend:            runtime switching between uniform and non-uniform modes
 *
 * Cost: always evaluates rau_warpf and performs the lerp, regardless of M.
 *   For fixed M known at compile time, use the NON_UNIFORM_VEL #define —
 *   the compiler eliminates the unused branch entirely.
 *   Use rau_sincos_mf only when M varies at runtime.
 *
 * Input:     RAU (0 to 4 = full revolution)
 * M:         blend factor ∈ [0.0, 1.0]
 * Max error: same as rau_sincosf at M=1; increases toward M=0 (no warp correction)
 */
void rau_sincos_mf(float input_t, float M, float *sin_out, float *cos_out);

/* ══════════════════════════════════════════════════════════════════════════
 * §2  INVERSE / PHASE-RECOVERY FUNCTIONS
 *
 * These recover RAU angle parameters from trigonometric ratios.
 * They form the inverse pipeline to rau_sincosf:
 * ══════════════════════════════════════════════════════════════════════════ */

/* rau_r_arcsinf — algebraic inverse from sin only.  Returns w ∈ [0,1].
 *
 * Solves  sin = w / sqrt((1−w)² + w²)  algebraically for w.
 * No polynomial approximation — uses only arithmetic and one sqrt.
 *
 * Singularity guard at 45° (sin = 1/√2, where sin = cos):
 *   The denominator 2sin²−1 passes through zero at this point.
 *   Guarded by returning 0.5 exactly — matches the warp fixed point warp(0.5)=0.5.
 *
 * err = 1: input outside [−1, 1] or non-finite.
 * Max error: 1.13e-14  (~46 effective bits — EXACT in practice)
 */
float rau_r_arcsinf(float s, int *err);

/* rau_r_arccosf — algebraic inverse from cos only.  Returns w ∈ [0,1].
 *
 * Solves  cos = (1−w) / sqrt((1−w)² + w²)  algebraically for w.
 * No polynomial approximation — uses only arithmetic and one sqrt.
 *
 * Use when only cos is available.
 * When both sin and cos are available, prefer rau_r_arctanf.
 *
 * Singularity guard at 45° (cos = 1/√2): returns 0.5 exactly.
 *
 * err = 1: input outside [−1, 1] or non-finite.
 * Max error: 1.13e-14  (~46 effective bits — EXACT in practice)
 */
float rau_r_arccosf(float c, int *err);

/* rau_invpolyf — fractional phase recovery.  w → fractional t ∈ [0,1).
 *
 * Recovers the fractional quadrant parameter t from the diagonal ratio w.
 * This is the approximate inverse of rau_warpf:
 *   given w = warp(t), returns t ≈ invpoly(w)
 *
 * Implementation uses rau_arctan_adjf with a reflection at w = 0.5:
 *   w ≥ 0.5:  t = (π/2 − arctan((1−w)/w)) · (2/π)   — complementary angle
 *   w < 0.5:  t =        arctan(w/(1−w))  · (2/π)
 * The reflection keeps the polynomial input always in [0, 0.5] where it is
 * well-conditioned, giving uniform accuracy across the full [0,1] domain.
 *
 * Note: this inverts warp(), not rau_sincosf().
 *   To recover a full angle φ ∈ [0,4) from (cos, sin), use rau_atan2f —
 *   it calls rau_r_arctanf and rau_invpolyf internally with quadrant reconstruction.
 *
 * err = 1: input outside [0, 1] or non-finite.
 * Max error: ~1.05e-4  (~13.2 effective bits, ~3 decimal places)
 */
float rau_invpolyf(float w, int *err);

/* rau_atan2f — full RAU phase recovery, analogous to atan2(y, x).
 * -----------------------------------------------------------------
 * Returns a full-circle phase φ in [0, 4), using the RAU convention:
 *   0   -> +x axis
 *   1   -> +y axis
 *   2   -> -x axis
 *   3   -> -y axis
 *
 * rau_atanf — principal-value inverse tangent, analogous to atan(u).
 * ------------------------------------------------------------------
 * Returns a signed value in [-1, 1], where:
 *   -1 -> -pi/2
 *    0 -> 0     [analogous to atan(x), normalized to RAU units]
 *   +1 -> +pi/2
 * instead of radians.
 *
 * Special cases:
 *   NaN      -> err = 1, return NAN
 *   +/-inf   -> +/-1
 *   +/-0     -> preserves signed zero
 */
float rau_atan2f(float ry, float rx, int *err);
float rau_atanf(float x, int *err);

/* Conversion helpers:
   signed radians = rau_atan2_signed_radians(rau_atan2f(y, x))
   signed degrees  = rau_atan2_signed_degrees(rau_atan2f(y, x))
   RAU to radians  = rau * (pi / 2)
   radians to RAU  = rad * (2 / pi)
 */
float rau_atan2_signed_radians(float phi_rau);
float rau_atan2_signed_degs(float phi_rau);

#ifdef __cplusplus
}
#endif

#endif /* RADICALTRIG_H */

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
