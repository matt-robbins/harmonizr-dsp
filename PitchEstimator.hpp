#ifndef _pitchest_
#define _pitchest_

#include <cstdio>
#include <string>
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

class PitchEstimator {
    public:
        virtual ~PitchEstimator() {}; // destructor, use it to call destructor of the inherit classes
        virtual float estimate(CircularAudioBuffer& b) = 0;
};

class PitchEstimatorYIN : public PitchEstimator {
public:
    PitchEstimatorYIN(int MaxT=800, int l2nfft=11, float threshold=0.5, int nmed=7);
    virtual ~PitchEstimatorYIN();

    virtual float estimate(CircularAudioBuffer& b);
    
private:
    int maxT;
    float threshold;

    int l2nfft = 11;
    int nfft = 0;
    int nmed = 7; // median filtering of results

    float *cmdf = nullptr;
    float *Tbuf = nullptr;
    float *Tsrt = nullptr;
    int Tix = 0;
    float fs = 44100;

    #ifdef __APPLE__
        FFTSetup fft_s;
        DSPSplitComplex fft_in, fft_in2, fft_out, fft_out2, fft_buf;
    #else
        kiss_fft_cfg fft_s, ifft_s;
        kiss_fft_cpx *fft_in, *fft_out, *fft_out2;
    #endif
};

class PitchEstimatorSHS : public PitchEstimator {
public:
    PitchEstimatorSHS(int MaxT=800, int l2nfft=11, float threshold=0.5, int nmed=7, float fs=44100.0f);
    virtual ~PitchEstimatorSHS();

    virtual float estimate(CircularAudioBuffer& b);
    
private:
    int maxT;
    float threshold;

    int n_oct;
    int points_per_oct = 48;
    float minf = 50; // this is Hz
    float maxf = 1600.0;
    float FS = 44100.0;

    int l2nfft = 11;
    int nfft = 0;
    int nspec = 0;
    int nmed = 7; // median filtering of results
    int npad = 2;

    float *linspec = nullptr;
    float *linspec_pad = nullptr;
    float *linpeaks_pad = nullptr;
    float *linpeaks = nullptr;
    float *logspec = nullptr;
    float *logidx = nullptr;
    float *Tbuf = nullptr;
    float *Tsrt = nullptr;
    int Tix = 0;
    Window w;

    #ifdef __APPLE__
        FFTSetup fft_s;
        DSPSplitComplex fft_in, fft_out, fft_buf;
    #else
        kiss_fft_cfg fft_s;
        kiss_fft_cpx *fft_in, *fft_out;
    #endif

};

#endif
