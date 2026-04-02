#include <cstdio>
#include <string>
#include <random>
#include "CircularAudioBuffer.hpp"
#include "Window.hpp"

#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include "kiss_fft.h"
#import <algorithm>
#include <android/log.h>

template <typename T>
T clamp(T input, T low, T high) {
    return std::min(std::max(input, low), high);
}

#endif

class PulseModel {
public:

    PulseModel(int l2nfft=11, float fc, float fs, int maxT);
    ~PulseModel() {
        free(fd_lpf);
    }

    void init(float fc, float fs);
    float * get_pulse_pointer();

    void get_minphase_pulse(float *data, int N, bool voiced, float T);
private:
    int l2nfft;
    int nfft;
    int maxT;

    int n_pulse{16}; // allows for 4 octaves of pitch increase

    #ifdef __APPLE__
        FFTSetup fft_s;
        DSPSplitComplex fft_in, fft_in2, fft_out, fft_out2, fft_buf, spec_env;
    #else
        kiss_fft_cfg fft_s, ifft_s;
        kiss_fft_cpx *fft_in, *fft_out, *fft_out2;
    #endif

    std::normal_distribution<float> dist{0.0f,1.0f};
    std::minstd_rand gen{0};
    Window w{Window::Hann,64};
    float * fd_lpf;

    float ** pulses;
    int pulse_ix{0};
};