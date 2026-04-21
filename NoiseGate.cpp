//
//  NoiseGate.cpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 4/1/26.
//
#include "NoiseGate.hpp"
#include <cmath>

NoiseGate::NoiseGate(float thresh_db, float hyst_db, float fs, float hold_s) :
thresh_db{thresh_db}, hyst_db{hyst_db}, fs{fs}, hold_s{hold_s} {
    thresh = std::pow(10.f, thresh_db / 20.f);
    hyst = std::pow(10.f, hyst_db / 20.f);
    hold_t = (int) hold_s*fs;
}

void NoiseGate::compute(float * out, float * in, int N) {
    for (int k = 0; k < N; k++){
        out[k] = compute_one(in[k]);
    }
}

void NoiseGate::set_thresh_db(float thresh_db) {
    this->thresh_db = thresh_db;
    thresh = std::pow(10.f, thresh_db / 20.f);
}

float NoiseGate::get_thresh_db() {
    return thresh_db;
}

float NoiseGate::get_gain() {
    return gain;
}

float NoiseGate::compute_one(float in) {
    env = env*0.999 + std::fabs(in)*0.001;
    //gain = gain*0.999 + target_gain*0.001;
    
    switch (state){
        case (Closed):
            //gain = gain*0.9999 + target_gain*0.0001;
            if (env >= thresh) {
                state = Open;
                target_gain = 1.0;
            }
            
        case (Open):
            //gain = gain*0.9999 + target_gain*0.0001;
            if (env >= thresh) {
                t = 0;
            }
            if (env < (thresh / hyst) && ++t > hold_t) {
                state = Closed;
                t = 0;
                target_gain = 0.0;
            }
    }
    //gain = gain*0.9999 + target_gain*0.0001;
    gain = gain_tracker.compute(target_gain);
        
    return in * gain;
}

float IirTracker::compute(float in) {
    
    float out = y[0]*2*p - y[1]*p*p + in*(1 - 2*p + p*p);
    
    y[1] = y[0]; y[0] = out;
    
    return out;
}
