//
//  SpectralProcessor.cpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 6/10/26.
//

#include "SpectralProcessor.hpp"

void SpectralProcessor::update_spectrum() {
    
    if (mFrozen) {
        // copy time-domain data into split-complex format for fft and Window!
        float * buf = b.getContiguousRelative(mNfft);

        // note: we can't use this cause we need to window!
        //vDSP_ctoz(reinterpret_cast<DSPComplex *>(buf), 2, &fft_x, 1, mNfft/2);

        for (int k = 0; k < mNfft/2; k++) {
            fft_x.realp[k] = buf[2*k] * w[2*k];
            fft_x.imagp[k] = buf[2*k+1] * w[2*k+1];
        }
        vDSP_fft_zrip(fft_s, &fft_x, 1, mP2nfft, FFT_FORWARD);
    }
    else {
        memset(fft_x.realp, 0, sizeof(float) * mNfft/2);
        memset(fft_x.imagp, 0, sizeof(float) * mNfft/2);
    }
    // store spectrum

    mag_spec[0] = fft_x.realp[0];
    mag_spec[mNfft/2] = fft_x.realp[1];
    
    for (int k = 1; k < mNfft/2; k++)
    {
        float r1,c1;
        r1 = fft_x.realp[k]; c1 = fft_x.imagp[k];
        float mag = sqrtf(r1*r1 + c1*c1) / (mNfft * (mLap-1));
        if (mFrozen && mWasFrozen)
            mag_spec[k] = mag_spec[k]; //*0.99 + 0.1*mag;
        else
            mag_spec[k] = mag;
        
        //mag_spec[k] =
        ph_spec[k] = atan2(r1, c1);
        // std::cout << "|X|=" << mag_spec[k] << "\n";
    }
    mWasFrozen = mFrozen;
}

void SpectralProcessor::generate() {
    
    // just set dc and nyquist to zero.
    fft_x.realp[0] = 0;
    fft_x.imagp[0] = 0;
    
    // keep magnitude spectrum, inject random phase, keep it conjugate symmetric
    for (int k = 1; k < mNfft/2; k++)
    {
        // random phase
        float phi = dist(rng);
        
        fft_x.realp[k] = mag_spec[k] * cosf(phi);
        fft_x.imagp[k] = mag_spec[k] * sinf(phi);
    }
    
    // inverse transform
    vDSP_fft_zrip(fft_s, &fft_x, 1, mP2nfft, FFT_INVERSE);
    
    // de-interleave and store in buffer
    vDSP_ztoc(&fft_x, 1, reinterpret_cast<DSPComplex *>(ob.get_current()), 2, mNfft / 2);
}

void SpectralProcessor::freeze(bool frozen) {
    mWasFrozen = mFrozen;
    mFrozen = frozen;
}

bool SpectralProcessor::frozen() {
    return mFrozen;
}

float SpectralProcessor::process_frame(float x) {
    b.pushValue(x);
    
    float y = ob.olapAdd();
    
    if (ob.check() >= 0) {
        update_spectrum();
        generate();
    }
    
    return y;
}
