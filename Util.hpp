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

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

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

//float cubic_interp(float *v, float ix)
//{
//    int iix = (int) floorf(ix);
//    return cubic(v + iix - 1, ix - iix);
//}
//
//float quadratic_peak(float *pv, float *v, int ix)
//{
//    float pix = 0.5*(v[ix-1] - v[ix+1])/(v[ix+1] + v[ix-1] - 2*v[ix]);
//    *pv = v[ix] - .25*(v[ix-1]-v[ix+1])*pix;
//    return pix;
//}
//
//float linear (float *v, float a)
//{
//    return v[0] * (1 - a) + v[1] * a;
//}
//
//float linear_interp(float *v, float ix)
//{
//    int iix = (int) floorf(ix);
//    return linear(v + iix, ix - iix);
//}

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
    
    if (p.realp == 0 || p.imagp == 0) {
        std::cout << "Memory allocation failure!\n";
    }
}

inline void fft_free(DSPSplitComplex &p) {
    free(p.realp);
    free(p.imagp);
}

class Counter {
public:
    Counter(int period) : m_period{period}, m_ix{0} {}
    
    bool update() {
        if (++m_ix >= m_period){
            m_ix = 0;
        }
        return m_ix == 0;
    }
    void setPeriod(int period) {
        m_period = period;
    }
private:
    int m_period = 10;
    int m_ix = 0;
};

template <typename T>
class StateMemVariable {
public:
    StateMemVariable(T init) : m_value{init} {}
    void set(T val) {
        m_last_value = m_value;
        m_value = val;
    }
    StateMemVariable& operator=(T other) {
        set(other);
        return *this;
    }
    
    operator T() const { return m_value; }
    auto operator<=>(T other) const { return m_value <=> other; }
    
    T value() {
        return m_value;
    }
    T prev() {
        return m_last_value;
    }
    T diff() {
        return m_value - m_last_value;
    }
    //constexpr auto operator<=>(const StateMemVariable&) const = default;
    
private:
    T m_value = (T) 0;
    T m_last_value = (T) 0;
};

// critically-damped smoothing operator. Pole position is p
class IirTracker {
public:
    IirTracker(float p=0.999f) : p{p} { };
    float compute(float in) {
        float out = y[0]*2*p - y[1]*p*p + in*(1 - 2*p + p*p);
        y[1] = y[0]; y[0] = out;
        return out;
    };
    void setP(float p) {
        this->p = p;
    }
    void setT(float t) {
        this->p = expf(-1.f/t);
    }
private:
    float y[2] = {0.f, 0.f};
    float p;
};

class RampTracker {
public:
    RampTracker(float speed=1.0) : m_speed{speed} {}
    float compute() {
        float diff = (m_target - m_current);
        if (fabs(diff) > m_speed)
            diff = sgn(diff) * m_speed;
        
        m_current += diff;
        return m_current;
    }
    void setTarget(float t) {
        m_target = t;
    }
private:
    float m_speed = 1.f;
    float m_target = 0.0;
    float m_current = 0.0;
};

inline float midi_note_to_T(float nn, float fs, float base=440.f, float tet=12.f) {
    float f = 440.f * powf(2.f, (nn - 69.f)/12.f);
    float T = fs/f;
    return T;
}

static constexpr float MIDI_A4 = 69.0;
inline float T_to_midi_note(float T, float Fs, float base=440.f, float tet=12.f) {
    return MIDI_A4 + log2f (Fs / (T * base)) * tet;
}

#endif
