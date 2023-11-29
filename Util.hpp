#ifndef util_h
#define util_h
#include <cmath>

inline float cubic (float *v, float a)
{
    float b, c;
    
    b = 1 - a;
    c = a * b;
    return (1.0f + 1.5f * c) * (v[1] * b + v[2] * a)
        - 0.5f * c * (v[0] * b + v[1] + v[2] + v[3] * a);
}

inline float cubic_a(float *v, float ix) {
    int ixi = floorf(ix);
    float f = ix - ixi;
    return cubic(v + ixi-1, f);
}

inline float cubic_v(float v0, float v1, float v2, float v3, float a) {
    float b, c;
    
    b = 1 - a;
    c = a * b;
    return (1.0f + 1.5f * c) * (v1 * b + v2 * a)
        - 0.5f * c * (v0 * b + v1 + v2 + v3 * a);
}

#endif
