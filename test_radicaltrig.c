#include <stdio.h>
#include <float.h>
#include "radicaltrig.h"
#include <math.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#ifndef M_PI_2
#define M_PI_2 1.57079632679489661923
#endif

static int failures = 0;

static void check(int cond, const char *name) {
    if (!cond) {
        ++failures;
        printf("FAIL: %s\n", name);
    }
}

static float fmax3(float a, float b, float c) {
    return fmaxf(a, fmaxf(b, c));
}

static int check_angle(float got, float expect, float tol, const char *name) {
    if (fabsf(got - expect) > tol) {
        printf("FAIL %s: got=%f expect=%f\n", name, got, expect);
        return 0;
    }
    return 1;
}

static int is_negzero(float x) {
    return x == 0.0f && signbit(x);
}

static int is_poszero(float x) {
    return x == 0.0f && !signbit(x);
}

static void test_atan2_special_cases(void)
{
    int e = 0;
    float a;

    struct {
        const char *name;
        float y;
        float x;
        float expect;
        int expect_nan;
    } cases[] = {
        { "(+inf,  0)",  INFINITY,   0.0f,  1.0f, 0 },
        { "(-inf,  0)", -INFINITY,   0.0f,  3.0f, 0 },
        { "(  0, +inf)", 0.0f,  INFINITY,   0.0f, 0 },
        { "(  0, -inf)", 0.0f, -INFINITY,   2.0f, 0 },
        { "(+inf,+inf)",  INFINITY,  INFINITY,  0.5f, 0 },
        { "(+inf,-inf)",  INFINITY, -INFINITY,  1.5f, 0 },
        { "(-inf,+inf)", -INFINITY,  INFINITY,  3.5f, 0 },
        { "(-inf,-inf)", -INFINITY, -INFINITY,  2.5f, 0 },
        { "NaN input",    NAN,       1.0f,      NAN,  1 },
        { "zero vector",  0.0f,      0.0f,      NAN,  1 },
    };

    for (size_t i = 0; i < sizeof(cases)/sizeof(cases[0]); ++i) {
        e = 0;
        a = rau_atan2f(cases[i].y, cases[i].x, &e);

        if (cases[i].expect_nan) {
            check(e != 0, cases[i].name);
            check(isnan(a), cases[i].name);
        } else {
            char buf[128];
            snprintf(buf, sizeof(buf), "%s err", cases[i].name);
            check(e == 0, buf);

            snprintf(buf, sizeof(buf), "%s angle", cases[i].name);
            check(fabsf(a - cases[i].expect) < 1e-7f, buf);
        }
    }
}

int main(void) {
    int e;
    float maxe;

    maxe = 0.0f;
    for (int i = 0; i <= 1000; ++i) {
        float t = 0.001f + (0.999f - 0.001f) * (float)i / 1000.0f;
        float w = rau_warpf(t);
        float ref = sinf(t * (float)M_PI_2) / (sinf(t * (float)M_PI_2) + cosf(t * (float)M_PI_2));
        maxe = fmaxf(maxe, fabsf(w - ref));
    }
    check(maxe < 5e-7f, "warpf");

    maxe = 0.0f;
    for (int i = 0; i <= 1000; ++i) {
        float p = 0.01f + (3.99f - 0.01f) * (float)i / 1000.0f;
        float s, c;
        rau_sincosf(p, &s, &c);
        float rs = sinf(p * (float)M_PI_2);
        float rc = cosf(p * (float)M_PI_2);
        maxe = fmax3(maxe, fabsf(s - rs), fabsf(c - rc));
    }
    check(maxe < 1e-6f, "sincosf");

    maxe = 0.0f;
    for (int i = 0; i <= 100; ++i) {
        float p = 0.05f + (0.95f - 0.05f) * (float)i / 100.0f;
        float t = rau_tanf(p);
        float ref = tanf(p * (float)M_PI_2);
        maxe = fmaxf(maxe, fabsf(t - ref));
    }
    check(maxe < 1e-4f, "tanf");

    maxe = 0.0f;
    for (int i = 0; i <= 500; ++i) {
        float p = 0.01f + (0.99f - 0.01f) * (float)i / 500.0f;
        float s, c;
        rau_sincosf(p, &s, &c);
        float w1 = rau_r_arctanf(s, c, &e);
        float wtrue = rau_warpf(p);
        maxe = fmaxf(maxe, fabsf(w1 - wtrue));
    }
    check(maxe < 1e-7f, "inverse w recovery");

    maxe = 0.0f;
    for (int i = 0; i <= 500; ++i) {
        float p = 0.01f + (0.49f - 0.01f) * (float)i / 500.0f;
        float s, c;
        rau_sincosf(p, &s, &c);
        float w2 = rau_r_arcsinf(s, &e);
        float w3 = rau_r_arccosf(c, &e);
        float wtrue = rau_warpf(p);
        maxe = fmaxf(maxe, fmaxf(fabsf(w2 - wtrue), fabsf(w3 - wtrue)));
    }
    check(maxe < 1e-5f, "principal inverse recovery");

    maxe = 0.0f;
    for (int i = 0; i <= 500; ++i) {
        float t = 0.01f + (0.99f - 0.01f) * (float)i / 500.0f;
        float w = rau_warpf(t);
        float tr = rau_invpolyf(w, &e);
        maxe = fmaxf(maxe, fabsf(tr - t));
    }
    check(maxe < 2e-4f, "invpoly");

    maxe = 0.0f;
    for (int i = 0; i <= 500; ++i) {
        float p = 0.05f + (3.95f - 0.05f) * (float)i / 500.0f;
        float s, c;
        rau_sincosf(p, &s, &c);
        float q = rau_atan2f(s, c, &e);
        float ref = atan2f(s, c) * (2.0f / (float)M_PI);
        if (ref < 0.0f) ref += 4.0f;
        float d = fabsf(q - ref);
        if (d > 2.0f) d = 4.0f - d;
        maxe = fmaxf(maxe, d);
    }
    check(maxe < 2e-4f, "full_phi");

    e = 0;
    (void)rau_r_arcsinf(2.0f, &e);
    check(e != 0, "domain error arcsin");

    e = 0;
    (void)rau_r_arccosf(-2.0f, &e);
    check(e != 0, "domain error arccos");

    e = 0;
    (void)rau_invpolyf(1.2f, &e);
    check(e != 0, "domain error invpoly");

    e = 0;
    (void)rau_atan2f(0.0f, 0.0f, &e);
    check(e != 0, "domain error full_phi");

    if (failures) {
        printf("%d test(s) failed\n", failures);
        return 1;
    }

    float a;

    a = rau_atan2f(0.0f, 0.0f, &e);
    //check(e != 0, "zero-vector error");
    check(isnan(a), "zero-vector nan");

    e = 0;
    a = rau_atan2f(NAN, 1.0f, &e);
    check(e != 0, "nan y error");

    e = 0;
    a = rau_atan2f(1.0f, NAN, &e);
    check(e != 0, "nan x error");

    e = 0;
    a = rau_atan2f(0.0f, 1.0f, &e);
    check(e == 0, "poszero east err");
    check(fabsf(a - 0.0f) < 1e-7f, "poszero east angle");

    e = 0;
    a = rau_atan2f(-0.0f, 1.0f, &e);
    check(e == 0, "negzero east err");
    check(is_poszero((a * (float)M_PI_2)), "negzero east signed-zero");

    e = 0;
    a = rau_atan2f(0.0f, -1.0f, &e);
    check(e == 0, "poszero west err");
    check(fabsf(a - 2.0f) < 1e-7f, "poszero west angle");

    e = 0;
    a = rau_atan2f(-0.0f, -1.0f, &e);
    check(e == 0, "negzero west err");
    check(fabsf(a - 2.0f) < 1e-7f, "negzero west angle");

    e = 0;
    a = rau_atan2f(INFINITY, 1.0f, &e);
    check(e == 0, "posinf north err");
    check(fabsf(a - 1.0f) < 1e-7f, "posinf north angle");

    e = 0;
    a = rau_atan2f(-INFINITY, 1.0f, &e);
    check(e == 0, "neginf south err");
    check(fabsf(a - 3.0f) < 1e-7f, "neginf south angle");

    e = 0;
    a = rau_atan2f(1.0f, INFINITY, &e);
    check(e == 0, "finite posinf x err");
    check(fabsf(a - 0.0f) < 1e-7f, "finite posinf x angle");

    e = 0;
    a = rau_atan2f(1.0f, -INFINITY, &e);
    check(e == 0, "finite neginf x err");
    check(fabsf(a - 2.0f) < 1e-7f, "finite neginf x angle");

    // Test for -inf and inf -- atan2 special cases
    test_atan2_special_cases();

    if (failures) {
        printf("%d test(s) failed\n", failures);
        return 1;
    }

    printf("All special-case tests passed\n");
    printf("All tests passed\n");

    int ok = 1;
    int err = 0;
    float phi;

    phi = rau_atan2f(0.0f, 1.0f, &err);
    ok &= check_angle(rau_atan2_signed_radians(phi), 0.0f, 1e-6f, "east signed");
    ok &= check_angle(rau_atan2_signed_degs(phi), 0.0f, 1e-5f, "east deg");

    phi = rau_atan2f(1.0f, 0.0f, &err);
    ok &= check_angle(rau_atan2_signed_radians(phi), (float)M_PI_2, 1e-6f, "north signed");
    ok &= check_angle(rau_atan2_signed_degs(phi), 90.0f, 1e-5f, "north deg");

    phi = rau_atan2f(0.0f, -1.0f, &err);
    ok &= check_angle(rau_atan2_signed_radians(phi), -(float)M_PI, 1e-6f, "west signed");
    ok &= check_angle(rau_atan2_signed_degs(phi), 180.0f, 1e-5f, "west deg");

    phi = rau_atan2f(-1.0f, 0.0f, &err);
    ok &= check_angle(rau_atan2_signed_radians(phi), -(float)M_PI_2, 1e-6f, "south signed");
    ok &= check_angle(rau_atan2_signed_degs(phi), 270.0f, 1e-5f, "south deg");

    return ok ? 0 : 1;
}
