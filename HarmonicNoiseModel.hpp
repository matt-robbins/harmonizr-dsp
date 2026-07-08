//
//  HarmonicNoiseModel.hpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 7/7/26.
//

#include <Accelerate/Accelerate.h>
#import <random>

#include "Util.hpp"
#include "Window.hpp"
#include "CircularAudioBuffer.hpp"

class HarmonicNoiseModel {
public:
    HarmonicNoiseModel(int p2nfft, int nH) : mP2nfft{p2nfft}, mNfft{0x01 << p2nfft}, b{p2nfft+1}, w{Window::Hann, 1<<p2nfft} {
        A.resize(nH);
        fft_s = vDSP_create_fftsetup(mP2nfft, FFT_RADIX2);

        fft_alloc(fft_x, mNfft);
        fft_alloc(A_model, max_harm);
        mag_spec = (float *) calloc(mNfft, sizeof(float));
    }
    
    void compute_model(float T, float * data, int N) {
        float f0 = (float)mNfft/T;
        fprintf(stderr, "f0 = %f\n", f0);
        int h_ix = 0;
        float ratio = 1.;
        int miss_cnt = 0;
        float h_max = 0.f;
        
        memset(fft_x.realp, 0, mNfft * sizeof(float));
        memset(fft_x.imagp, 0, mNfft * sizeof(float));
        
        for (int k = 0; k < 3*T; k+=2)
        {
            fft_x.realp[k] = data[k] * w.value((float) k / 3*T);
            fft_x.imagp[k] = data[k+1] * w.value((float) (k+1) / 3*T);
        }
        
        vDSP_fft_zrip(fft_s, &fft_x, 1, mP2nfft, FFT_FORWARD);
        vDSP_zvmags(&fft_x, 1, mag_spec, 1, mNfft/2);
        
        float fix = 0;
        
        // iterate up through harmonics, or until Fs/16
        while (fix < mNfft/8 && h_ix < max_harm)
        {
            fix += f0;
            h_ix++;
            //fprintf(stderr, "fix=%f\n",fix);
            float max = 0.f;
            int max_ix = 0;
            
            // find the max magnitude around the next likely harmonic
            for (int k = floorf(fix - f0/4); k < ceilf(fix + f0/4); k++)
            {
                if (mag_spec[k] > max)
                {
                    max = mag_spec[k]; max_ix = k;
                }
            }
            // if it's really a peak...
            if (max > mag_spec[max_ix+1] && max > mag_spec[max_ix-1])
            {
                float pv;
                float pix = max_ix + quadratic_peak(&pv, mag_spec, max_ix);
                
                if (pv > h_max) h_max = pv;
                
                //fprintf(stderr, "px = %f\n",pix);
                //fprintf(stderr, "pv = %f\n",pv);
                
                f0 = f0 * (1-ratio) + ratio * (float)(pix)/h_ix;
                A_model.realp[h_ix] = linear_interp(fft_x.realp, pix);
                A_model.imagp[h_ix] = linear_interp(fft_x.imagp, pix);
                fprintf(stderr, "A[%d] = %f\n", h_ix, sqrtf(A_model.realp[h_ix]*A_model.realp[h_ix] + A_model.imagp[h_ix]*A_model.imagp[h_ix]));
            }
            else {
                if(++miss_cnt > 3)
                    break;
            }
            
            ratio /= 1.2;
        }
        
        fprintf(stderr,"f0_hat = %f\n\n", f0);
    }
private:
    std::vector<float> A;
    Window w;
    CircularAudioBuffer b;
    
    FFTSetup fft_s;
    DSPSplitComplex fft_x;
    DSPSplitComplex A_model;
    
    float * mag_spec;
    
    int mP2nfft;
    int mNfft;
    int max_harm = 30;
};
