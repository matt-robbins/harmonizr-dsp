//
//  MinphaseModel.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 7/7/26.
//

#include <Accelerate/Accelerate.h>
#include <random>
#include "Util.hpp"
#include "Window.hpp"

class MinphaseModel {
public:
    MinphaseModel(int p2nfft, int olap, float sampleRate) : mP2nfft{p2nfft}, mNfft{0x01 << p2nfft}, w{Window::Hann, 1 << p2nfft} {
        fft_s = vDSP_create_fftsetup(mP2nfft, FFT_RADIX2);
        fft_alloc(fft_in, mNfft);
        fft_alloc(fft_out, mNfft);
        fft_alloc(fft_buf, mNfft);
        
        synth_pulse = (float **) calloc(n_synth_pulse, sizeof(float *));
        
        for (int k = 0; k < n_synth_pulse; k++)
            synth_pulse[k] = (float *) calloc(2*maxT, sizeof(float));
        
        fd_lpf = (float *) calloc(mNfft, sizeof(float));
        
        int fc = (int) ((4000. / sampleRate) * mNfft); // cutoff frequency chosen to preserve formants
        int bw = fc / 2;
        //fprintf(stderr, "fc = %d\n", fc);
        for (int k = 1; k < mNfft/2; k++)
        {
            if (k <= fc-bw)
                fd_lpf[k] = fd_lpf[mNfft-k] = 1.0;
            else if (k < fc+bw)
                fd_lpf[k] = fd_lpf[mNfft-k] = 0.5 - 0.5 * cos(M_PI * (k - (fc - bw)) / (2*bw + 1));
            //fprintf(stderr, "%f\n", fd_lpf[k]);
        }

    }
    
    void calculate(float * input, int N, float T, bool voiced) {
        memset(fft_in.realp, 0, mNfft * sizeof(float));
        memset(fft_in.imagp, 0, mNfft * sizeof(float));
        static int count = 0;
        for (int k = 0; k < mNfft; k++)
        {
            //int ix = (start_ix + k) & cmask;
            fft_in.realp[k] = input[k] * w.value((float)k/(mNfft));
        }
        
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, mP2nfft, 1);
        // log(abs(fft))
        for (int k = 0; k < mNfft; k++)
        {
            float r1 = fft_out.realp[k];
            float c1 = fft_out.imagp[k];
            float mag = (r1*r1 + c1*c1);
            
            fft_in.realp[k] = logf(mag+0.00001)/2; // factor of 2 is square root
            fft_in.imagp[k] = 0;
        }

        // cepstrum
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, mP2nfft, -1);
        
        // fold and window
        fft_in.realp[0] = fft_out.realp[0]/mNfft;

        int cutoff = 300;
        if (voiced)
        {
            cutoff = roundf(T * 0.7);
        }
        for (int k = 1; k < cutoff; k++)
        {
            fft_in.realp[k] = 0.f;
            fft_in.imagp[k] = 0;
        }
        for (int k = cutoff; k < mNfft; k++)
        {
            fft_in.realp[k] = (fft_out.realp[k]*2)/mNfft;
            fft_in.imagp[k] = 0;
        }
        
        // fft_out contains smooth envelope (fft)
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, mP2nfft, 1);
        
        for (int k = 0; k < mNfft; k++)
        {
            fft_in.realp[k] = dist(rng); fft_in.imagp[k] = 0;
        }
        
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, mP2nfft, 1);
        // fft_out2 contains random spectrum, convolve with envelope
        for (int k = 0; k < mNfft; k++)
        {
            float ex = expf(fft_out.realp[k]);
            float flt = 1; //(1 + 9*(1 - fd_lpf[k]))/10;
            if (voiced)
            {
                flt = 1 - fd_lpf[k];
            }
            fft_in.realp[k] = flt * ex * fft_out2.realp[k]/mNfft;
            fft_in.imagp[k] = flt * ex * fft_out2.imagp[k]/mNfft;
        }
        
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, mP2nfft, -1);
        // fft_out2 now contains noise matching original sound spectrum
        
        for (int k = 0; k < mNfft; k++)
        {
            spec_env.realp[k] = 0.9*spec_env.realp[k] + 0.1*fft_out.realp[k];
            spec_env.imagp[k] = 0.9*spec_env.imagp[k] + 0.1*fft_out.realp[k];

            float ex = expf(spec_env.realp[k]);
            //float rnd = nd(de);
            
            if (voiced)
            {
                ex *= fd_lpf[k];
            }
            
            float fr = fabs(1. - ((float) k / (float) (1 + mNfft/2)));

            if (!voiced)
            {
                fr = 0.;
            }
            
            fft_in.realp[k] = ex * cos(spec_env.realp[k]);
            fft_in.imagp[k] = ex * sin(spec_env.realp[k]);
        }
        
        // get impulse response (ifft)
        vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, mP2nfft, -1); //inverse
        
        if (++synth_pulse_ix >= n_synth_pulse)
            synth_pulse_ix = 0;
        
        if (voiced && count < 10)
        {
            count++;
        }
        
        if (!voiced && count > 0)
        {
            count--;
        }
        
        float mix = (float) count / 100.0;
        //fprintf(stderr, "mix = %f\n", mix);
        
        for (int k = 0; k < 2*maxT; k++)
        {
            synth_pulse[synth_pulse_ix][k] = mix*fft_out.realp[k]/sqrtf(mNfft)/16;
            //synth_pulse[synth_pulse_ix][k] += fft_out2.realp[k]/16;
        }
    }
private:
    int mP2nfft;
    int mNfft;
    
    FFTSetup fft_s;
    DSPSplitComplex fft_in, fft_out, fft_out2, fft_buf, A_model, Hann, spec_env;
    Window w;
    float * mag_spec;
    float * ph_spec;
    float * fd_lpf;
    float ** synth_pulse;
    int n_synth_pulse = 10;
    int synth_pulse_ix = 0;
    int maxT = 800;
    float ** out;
    std::mt19937 rng;
    std::normal_distribution<float> dist{0, 1};
    
};
