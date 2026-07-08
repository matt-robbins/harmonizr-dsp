//
//  SpectralProcessor.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 6/10/26.
//

#include <Accelerate/Accelerate.h>
#import <random>

#include "Window.hpp"
#include "CircularAudioBuffer.hpp"
#include "Util.hpp"

class OlapAddBuffer {
public:
    OlapAddBuffer(int N = 2, int L = 256) : mN{N}, mL{L}, w{Window::Hann, L}, buf(N,std::vector<float>(L,0)), ix(N,0) {
        for (int k = 0; k < mN; k++) {
            ix[k] = (k * mL)/mN;
            std::cout << "ix[" << k << "] = " << ix[k] << "\n";
        }
    }
    ~OlapAddBuffer() {
    }
    
    void write(float * src, int L) {
        if (L != mL) {
            std::cout << "error!\n";
        }
        memcpy(buf[bufix].data(), src, L * sizeof(float));
    }
    
    float * get_current() {
        return buf[bufix].data();
    }

    float olapAdd() {
        float o = 0.f;
        for (int k = 0; k < mN; k++) {
            o += buf[k][ix[k]] * w[ix[k]];
        }
        increment();
        return o;
    }
    int check() {
        return next;
    }
private:
    Window w;
    std::vector<std::vector<float>> buf;
    std::vector<int> ix;
    
    int bufix = 0;
    int next = -1;
    int mN;
    int mL;
    
    void increment() {
        next = -1;
        for (int k = 0; k < mN; k++) {
            if (++ix[k] >= mL) {
                ix[k] = 0;
                next = k;
                bufix = k;
            }
        }
    }
};

class SpectralProcessor {
public:
    SpectralProcessor(int p2nfft, int olap) : mP2nfft{p2nfft}, mNfft{0x01 << p2nfft}, b{p2nfft+1}, ob{olap, 0x01 << (p2nfft)}, w{Window::Hann,1<<p2nfft}, rng{std::random_device{}()}, dist{0.0, M_PI*2.f}
    {
        fft_s = vDSP_create_fftsetup(mP2nfft, FFT_RADIX2);

        fft_alloc(fft_x, mNfft);
        
        mag_spec = (float *) calloc(mNfft, sizeof(float));
        ph_spec = (float *) calloc(mNfft, sizeof(float));

        mLap = olap;
    }
    SpectralProcessor(const SpectralProcessor& other) : mP2nfft{other.mP2nfft}, mNfft{other.mNfft}, b{other.mP2nfft+1}, ob{other.mLap, other.mNfft}, w{Window::Hann, other.mNfft}, rng{std::random_device{}()}, dist{0.0, M_PI*2.f} {
        fft_s = vDSP_create_fftsetup(mP2nfft, FFT_RADIX2);

        fft_alloc(fft_x, mNfft);
        
        mag_spec = (float *) calloc(mNfft, sizeof(float));
        ph_spec = (float *) calloc(mNfft, sizeof(float));
    }
    ~SpectralProcessor() {
        
        fft_free(fft_x);
        //fft_free(fft_X);
        
        free(mag_spec);
        free(ph_spec);
        
        vDSP_destroy_fftsetup(fft_s);
        
    }; // destructor, use it to call destructor of the inherit classes
    float process_frame(float x);
    void update_spectrum();
    
    void generate();
    void freeze(bool frozen);
    bool frozen();

private:
    OlapAddBuffer ob;
    CircularAudioBuffer b;
    Window w;
    
    int mP2nfft = 11;
    int mNfft = 0;
    int mLap = 0;
    
    bool mFrozen = false;
    bool mWasFrozen = false;
    
    FFTSetup fft_s;
    DSPSplitComplex fft_x;
    float * mag_spec;
    float * ph_spec;
    float ** out;
    
    std::mt19937 rng;
    std::uniform_real_distribution<float> dist;
};
