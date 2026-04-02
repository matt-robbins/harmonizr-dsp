//
//  NoiseGate.cpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 4/1/26.
//
#include "NoiseGate.hpp"
#include <cmath>

NoiseGate::NoiseGate(float thresh_db, float hyst_db, int fs, int hold_t) :
thresh_db{thresh_db}, hyst_db{hyst_db}, fs{fs}, hold_t{hold_t} {
    thresh = std::pow(10.f, thresh_db / 20.f);
    hyst = std::pow(10.f, hyst_db / 20.f);
}

void NoiseGate::compute(float * out, float * in, int N) {
    for (int k = 0; k < N; k++){
        out[k] = compute_one(in[k]);
    }
}

float NoiseGate::get_gain() {
    return gain;
}

float NoiseGate::compute_one(float in) {
    env = env*0.999 + std::fabs(in)*0.001;
    //gain = gain*0.999 + target_gain*0.001;
    
    switch (state){
        case (Closed):
            gain = gain*0.9999 + target_gain*0.0001;
            if (env >= thresh) {
                state = Open;
                target_gain = 1.0;
            }
            
        case (Open):
            gain = gain*0.999 + target_gain*0.001;
            if (env >= thresh) {
                t = 0;
            }
            if (env < (thresh / hyst) && ++t > hold_t) {
                state = Closed;
                t = 0;
                target_gain = 0.0;
            }
    }
        
    return in * gain;
}
