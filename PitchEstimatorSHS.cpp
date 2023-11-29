#include "PitchEstimator.hpp"
#include "Util.hpp"
#include <iostream>

static void fft_alloc(DSPSplitComplex &p, int nfft) {
    p.realp = (float *) calloc(nfft, sizeof(float));
    p.imagp = (float *) calloc(nfft, sizeof(float));
}

static void fft_free(DSPSplitComplex &p) {
    free(p.realp);
    free(p.imagp);
}

PitchEstimatorSHS::PitchEstimatorSHS(int MaxT, int l2nfft, float thresh, int nmed, float fs) :
maxT{MaxT}, l2nfft{l2nfft}, threshold{thresh}, nmed{nmed}, FS{fs},
w {Window(Window::Hamm, 0x01 << l2nfft)} {

    nfft = 0x01 << l2nfft;
    nspec = nfft/2 + 1;
    
    fft_s = vDSP_create_fftsetup(l2nfft, 2);

    fft_alloc(fft_in, nfft);
    fft_alloc(fft_out, nfft);
    fft_alloc(fft_buf, nfft);

    linspec_pad = (float *) calloc(nfft+npad*2, sizeof(float));
    linspec = linspec_pad + npad;
    linpeaks_pad = (float *) calloc(nfft+npad*2, sizeof(float));
    linpeaks = linpeaks_pad + npad;
    
    n_oct = ceilf( log2(maxf/minf));
    
    logspec = (float *) calloc(n_oct*points_per_oct, sizeof(float));
    logidx = (float *) calloc(n_oct*points_per_oct, sizeof(float));

    {
        int idx;
        float f;
        for (f = minf, idx = 0; idx < points_per_oct*n_oct; f*=powf(2,1.f/points_per_oct), idx++) {
            logidx[idx] = f*nfft/FS;
        }
    }
}

PitchEstimatorSHS::~PitchEstimatorSHS(){
    vDSP_destroy_fftsetup(fft_s);
    fft_free(fft_in);
    fft_free(fft_out);
    fft_free(fft_buf);
    free(linspec_pad);
    free(linpeaks);
    free(logspec);
    free(Tbuf);
    free(Tsrt);
}

float PitchEstimatorSHS::estimate(CircularAudioBuffer &b) {
    memset(fft_in.realp, 0, nfft * sizeof(float));
    memset(fft_in.imagp, 0, nfft * sizeof(float));
    
    b.copyRange(nfft, nfft, fft_in.realp);
    for (int i = 0; i < nfft; i++){
        fft_in.realp[i] *= w[i];
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);

    // compute amplitude
    for (int i = 0; i < nspec - 1; i++){
        linspec[i] = sqrtf(fft_out.realp[i]*fft_out.realp[i] + fft_out.imagp[i]*fft_out.imagp[i]);
    }
    // mask for points that are closer than 2 indices away from a local max
    memset(linpeaks, 0, nfft*sizeof(float));
    for (int i = 1; i < nfft - 1; i++){
        if (linspec[i] > linspec[i-1] && linspec[i] > linspec[i+1]){
            for (int j = -2; j <= 2; j++){
                linpeaks[i+j] = 1.0;
            }
        }
    }
    for (int i = 1; i < nspec - 1; i++){
        linpeaks[i] = linspec[i]*linpeaks[i];
    }
    // smoothing filter
    for (int i = 0; i < nspec; i++) {
        linspec[i] = 0.25*linpeaks[i-1] + 0.5*linpeaks[i] + 0.25*linpeaks[i+1];
    }

    // log interpolation
    for (int i = 0; i < n_oct * points_per_oct; i++) {
        logspec[i] = cubic_a(linspec, logidx[i]);
    }

    return 0.0f;
}
