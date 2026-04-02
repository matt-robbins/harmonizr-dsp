#ifndef util_h
#define util_h
#include <cmath>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <vector>

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#endif

#ifdef DEBUG
#  define D(x) x
#else
#  define D(x)
#endif // DEBUG

inline float cubic (float *v, float a)
{
    float b, c;
    
    b = 1 - a;
    c = a * b;
    return (1.0f + 1.5f * c) * (v[1] * b + v[2] * a)
        - 0.5f * c * (v[0] * b + v[1] + v[2] + v[3] * a);
}

inline float cubic(std::vector<float>::iterator v, float a) {
    float b, c;
    
    b = 1 - a;
    c = a * b;
    return (1.0f + 1.5f * c) * (v[1] * b + v[2] * a)
        - 0.5f * c * (v[0] * b + v[1] + v[2] + v[3] * a);
}

// inline float cubic (std::span<float> v, float a) {
//     float b, c;
    
//     b = 1 - a;
//     c = a * b;
//     return (1.0f + 1.5f * c) * (v[1] * b + v[2] * a)
//         - 0.5f * c * (v[0] * b + v[1] + v[2] + v[3] * a);
// }


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

// this is SUPER Dangerous -- be sure you know WTF you're doing
inline float valueAtIndexInterp(float *data, float ix) {
    float frac = ix - floorf(ix);
    int i = (int) floorf(ix-1);
    float u = cubic (data + i, frac);
    return u;
}

inline float median_idx(float *data, int n) {
    float min = *std::min_element(data,data+n);
    float sum = std::accumulate(data,data+n,-min*n);

    float csum = 0, last_csum = 0;
    int med_ix = 0;
    float median = 0;
    for (int k = 0; k < n; k++){
        last_csum = csum;
        csum += (data[k] - min);
        if (csum > sum/2){
            med_ix = k;
            break;
        }
    }

    if (sum == 0) {
        D(std::cout << "sum==0!" << std::endl;)
        return (float)n/2;
    }
    // linearly interpolate
    float frac = ((sum/2)-last_csum)/(csum-last_csum);
    median = med_ix+frac-0.5;
    return median;
}

inline void fft_alloc(DSPSplitComplex &p, int nfft) {
    p.realp = (float *) calloc(nfft, sizeof(float));
    p.imagp = (float *) calloc(nfft, sizeof(float));
}

inline void fft_free(DSPSplitComplex &p) {
    free(p.realp);
    free(p.imagp);
}

#endif
