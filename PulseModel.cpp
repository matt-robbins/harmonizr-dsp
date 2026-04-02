#include "PulseModel.hpp"
#include "Util.hpp"
#include <random>

PulseModel::PulseModel(int l2nfft, float fc, float fs, int maxT) 
    : l2nfft(l2nfft), maxT{maxT} {
    nfft = 2 << l2nfft;
    fd_lpf = (float *) calloc(nfft, sizeof(float));

    fft_s = vDSP_create_fftsetup(l2nfft, 2);
    fft_alloc(fft_in, nfft);
    fft_alloc(fft_out, nfft);
    fft_alloc(spec_env, nfft);

    pulses = (float **) calloc(n_pulse, sizeof(float*));

    for (int i = 0; i < n_pulse; i++) {
        pulses[i] = (float *) calloc(2*maxT, sizeof(float));
    }
    
    std::default_random_engine de(0); //seed
    init(fc,fs);
}

PulseModel::~PulseModel() {
    free(fd_lpf);
    fft_free(fft_in);
    fft_free(fft_out);
    fft_free(spec_env);
    vDSP_destroy_fftsetup(fft_s);

    for (int i = 0; i < n_pulse; i++) {
        free(pulses[i]);
    }
    free(pulses);
}

void PulseModel::init(float fc, float fs) {
    int fc = (int) ((fc / fs) * nfft); // cutoff frequency chosen to preserve formants
    int bw = fc / 2;
    //fprintf(stderr, "fc = %d\n", fc);
    for (int k = 1; k < nfft/2; k++)
    {
        if (k <= fc-bw)
            fd_lpf[k] = fd_lpf[nfft-k] = 1.0;
        else if (k < fc+bw)
            fd_lpf[k] = fd_lpf[nfft-k] = 0.5 - 0.5 * cos(M_PI * (k - (fc - bw)) / (2*bw + 1));
        //fprintf(stderr, "%f\n", fd_lpf[k]);
    }
}

void PulseModel::get_minphase_pulse(float *data, int N, bool voiced, float T){
    memset(fft_in.realp, 0, nfft * sizeof(float));
    memset(fft_in.imagp, 0, nfft * sizeof(float));
    static int count = 0;

    if (N > nfft) N = nfft;
    
    for (int k = 0; k < N; k++){
        fft_in.realp[k] = data[k] * w.value((float)k/(N));
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, l2nfft, 1);
    // log(abs(fft))
    for (int k = 0; k < nfft; k++)
    {
        float r1,c1;
        r1 = fft_out.realp[k]; c1 = fft_out.imagp[k];
        float mag = (r1*r1 + c1*c1);
        
        fft_in.realp[k] = logf(mag+0.00001)/2; // factor of 2 is square root
        fft_in.imagp[k] = 0;
    }

    // cepstrum
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, l2nfft, -1);
    
    // fold and window
    fft_in.realp[0] = fft_out.realp[0]/nfft;

    int cutoff = 300; // time domain samples
    if (voiced)
    {
        cutoff = roundf(T * 0.7);
    }
    for (int k = 1; k < cutoff; k++)
    {
        fft_in.realp[k] = 0.f;
        fft_in.imagp[k] = 0.f;
    }
    for (int k = cutoff; k < nfft; k++)
    {
        fft_in.realp[k] = (fft_out.realp[k]*2)/nfft;
        fft_in.imagp[k] = 0.f;
    }
    
    // fft_out contains smooth envelope (fft)
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, 11, 1);
    
    // generate random time-domain data
    for (int k = 0; k < nfft; k++)
    {
        fft_in.realp[k] = dist(gen); fft_in.imagp[k] = 0;
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, 1);
    // fft_out2 contains random spectrum, convolve with envelope
    for (int k = 0; k < nfft; k++)
    {
        float ex = expf(fft_out.realp[k]);
        float flt = 1; //(1 + 9*(1 - fd_lpf[k]))/10;
        if (voiced){
            flt = 1 - fd_lpf[k];
        }
        fft_in.realp[k] = flt * ex * fft_out2.realp[k]/nfft;
        fft_in.imagp[k] = flt * ex * fft_out2.imagp[k]/nfft;
    }
    
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out2, 1, &fft_buf, 11, -1);
    // fft_out2 now contains noise matching original sound spectrum
    
    // do time smoothing of spectral envelope and raise to e(th) power
    // e^(re + im*j) = e^re*e^(im*j) = e^re*cos(im) + j*e^re*sin(im))
    // re = e^re*cos(im)
    // im = e^re*sin(im)
    // Filter in the frequency domain by multiplying by our fd_lpf

    float smooth = 0.1;
    float smooth_c = 1 - smooth;
    for (int k = 0; k < nfft; k++)
    {
        spec_env.realp[k] = smooth_c*spec_env.realp[k] + smooth*fft_out.realp[k];
        spec_env.imagp[k] = smooth_c*spec_env.imagp[k] + smooth*fft_out.imagp[k];

        float ex = expf(spec_env.realp[k]);
        
        if (voiced){
            ex *= fd_lpf[k];
        }
        
        fft_in.realp[k] = ex * cos(spec_env.imagp[k]);
        fft_in.imagp[k] = ex * sin(spec_env.imagp[k]);
    }
    
    // get impulse response (ifft)
    vDSP_fft_zopt(fft_s, &fft_in, 1, &fft_out, 1, &fft_buf, l2nfft, -1); //inverse
    
    if (++pulse_ix >= n_pulse)
        pulse_ix = 0;
    
    if (voiced && count < 10){
        count++;
    }
    
    if (!voiced && count > 0){
        count--;
    }
    
    float mix = (float) count / 100.0;
    mix = 1.0f;
    //fprintf(stderr, "mix = %f\n", mix);
    
    for (int k = 0; k < 2*maxT; k++)
    {
        pulses[pulse_ix][k] = mix*fft_out.realp[k]/sqrtf(nfft);
        //pulses[pulse_ix][k] += fft_out2.realp[k]/16;
    }
            
}