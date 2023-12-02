#include "PitchEstimator.hpp"
#include <iostream>

static void fft_alloc(DSPSplitComplex &p, int nfft) {
    p.realp = (float *) calloc(nfft, sizeof(float));
    p.imagp = (float *) calloc(nfft, sizeof(float));
}

static void fft_free(DSPSplitComplex &p) {
    free(p.realp);
    free(p.imagp);
}

PitchEstimatorYIN::PitchEstimatorYIN(int MaxT, int l2nfft, float thresh, int nmed) : 
    maxT{MaxT}, l2nfft{l2nfft}, threshold{thresh}, nmed{nmed} {

    nfft = 0x01 << l2nfft;
    fft_s = vDSP_create_fftsetup(l2nfft, 2);

    fft_alloc(fft_in, nfft);
    fft_alloc(fft_in2, nfft);
    fft_alloc(fft_out, nfft);
    fft_alloc(fft_out2, nfft);
    fft_alloc(fft_buf, nfft);

    cmdf = (float *) calloc(maxT, sizeof(float));
    Tbuf = (float *) calloc(nmed, sizeof(float));
    Tsrt = (float *) calloc(nmed, sizeof(float));
}

PitchEstimatorYIN::~PitchEstimatorYIN(){
    vDSP_destroy_fftsetup(fft_s);
    fft_free(fft_in);
    fft_free(fft_in2);
    fft_free(fft_out);
    fft_free(fft_out2);
    fft_free(fft_buf);
    free(cmdf);
    free(Tbuf);
    free(Tsrt);
}

float PitchEstimatorYIN::estimate(CircularAudioBuffer &b) {

    static float old_period = 0;
    static int dead_count = 0;
    static int live_count = 0;
    
    //fprintf(stderr, "estimate pitch : imagp = %p\n", fft_in.imagp);
    memset(fft_in.realp, 0, nfft * sizeof(float));
    memset(fft_in.imagp, 0, nfft * sizeof(float));
    
    // copy first half, and compute FFT
    b.copyRange(2*maxT,maxT, fft_in.realp);
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);

    // copy the rest and compute again for large window
    b.copyRange(maxT,maxT, fft_in.realp + maxT);
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, 1);
    
    // conjugate small window and correlate with large window
    for (int k = 0; k < nfft; k++)
    {
        float r1,c1,r2,c2;
        r1 = fft_out.realp[k]; c1 = -fft_out.imagp[k];
        r2 = fft_out2.realp[k]; c2 = fft_out2.imagp[k];
        
        fft_in2.realp[k] = (r1*r2 - c1*c2);
        fft_in2.imagp[k] = (r1*c2 + r2*c1);
    }
    // inverse transform
    vDSP_fft_zopt(fft_s, &fft_in2, 1, &fft_out, 1, &fft_buf, 11, -1);
    
    float sumsq_ = fft_out.realp[0]/nfft;
    float sumsq = sumsq_;
    
    float df,cmdf1,cmdf2, sum = 0;
    
    float period = 0.0;
    
    // compute cumulative mean difference function
    cmdf2 = cmdf1 = 1;
    cmdf[0] = 1;
    for (int k = 1; k < maxT; k++)
    {
        float v1 = fft_in.realp[k]; // fft_in.realp still holds the raw time series data
        float v2 = fft_in.realp[k+maxT];
        
        sumsq -= v1*v1;
        sumsq += v2*v2;
        
        df = sumsq + sumsq_ - 2 * fft_out.realp[k]/nfft;
        sum += df;
        cmdf2 = cmdf1; cmdf1 = cmdf[k-1];
        cmdf[k] = (df * k) / sum;

        if (k > 0 && cmdf2 > cmdf1 && cmdf1 < cmdf[k] && cmdf1 < threshold && k > 20)
        {
            dead_count = 0;
            period = (float) (k-1) + 0.5*(cmdf2 - cmdf[k])/(cmdf2 + cmdf[k] - 2*cmdf1);
            
            break;
        }
    }
    
    if (period == 0 && dead_count < 10){
        live_count = 0;
        dead_count++;
        period = old_period;
    }
    else{
        live_count++;
    }
    
    old_period = period;
    Tbuf[Tix++] = period;
    
    //fprintf(stderr, "%f\n", sumsq);
    
    if (Tix >= nmed)
        Tix = 0;
    
    memcpy(Tsrt, Tbuf, nmed * sizeof(float));
    vDSP_vsort(Tsrt, (vDSP_Length) nmed, 1);
    
    return Tsrt[nmed/2];
}

// KISS-FFT port for Android

#if 0
float HarmonizerDSPKernel::estimate_pitch(int start_ix) {
    memset(fft_in, 0, nfft * sizeof(kiss_fft_cpx));

    for (int k = 0; k < maxT; k++) {
        int ix = (start_ix + k) & cmask;
        fft_in[k].r = cbuf[ix];
    }

    kiss_fft(fft_s, fft_in, fft_out);

    //memset(fft_in, 0, nfft * sizeof(kiss_fft_cpx));

    for (int k = maxT; k < 2 * maxT; k++) {
        int ix = (start_ix + k) & cmask;
        fft_in[k].r = cbuf[ix];
    }

    kiss_fft(fft_s, fft_in, fft_out2);

    // conjugate small window and correlate with large window
    for (int k = 0; k < nfft; k++) {
        float r1, c1, r2, c2;
        r1 = fft_out[k].r;
        c1 = -fft_out[k].i;
        r2 = fft_out2[k].r;
        c2 = fft_out2[k].i;

        fft_in[k].r = fd_lpf[k] * (r1 * r2 - c1 * c2);
        fft_in[k].i = fd_lpf[k] * (r1 * c2 + r2 * c1);
    }
    // inverse transform
    kiss_fft(ifft_s, fft_in, fft_out);

    float sumsq_ = fft_out[0].r / nfft;
    float sumsq = sumsq_;

    float df, cmdf, cmdf1, cmdf2, sum = 0;

    float period = 0.0;

    cmdf2 = cmdf1 = cmdf = 1;
    for (int k = 1; k < maxT; k++) {
        int ix1 = (start_ix + k) & cmask;
        int ix2 = (start_ix + k + maxT) & cmask;

        sumsq -= cbuf[ix1] * cbuf[ix1];
        sumsq += cbuf[ix2] * cbuf[ix2];

        df = sumsq + sumsq_ - 2 * fft_out[k].r / nfft;
        sum += df;
        cmdf2 = cmdf1;
        cmdf1 = cmdf;
        cmdf = (df * k) / sum;

        if (k > 0 && cmdf2 > cmdf1 && cmdf1 < cmdf && cmdf1 < threshold && k > 20) {
            period = (float) (k - 1) + 0.5 * (cmdf2 - cmdf) / (cmdf2 + cmdf - 2 * cmdf1);
            break;
        }
    }

    Tbuf[Tix++] = period;

    if (Tix >= nmed)
        Tix = 0;

    memcpy(Tsrt, Tbuf, nmed * sizeof(float));
    std::sort(Tsrt, Tsrt+nmed);
    return Tsrt[nmed / 2];
}
#endif
