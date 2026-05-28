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

int WinMain(int argc, char *argv[]) {
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
        float q = rau_full_phif(s, c, &e);
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
    (void)rau_full_phif(0.0f, 0.0f, &e);
    check(e != 0, "domain error full_phi");

    if (failures) {
        printf("%d test(s) failed\n", failures);
        return 1;
    }

    printf("All tests passed\n");
    return 0;
}
