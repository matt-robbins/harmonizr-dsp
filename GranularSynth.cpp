#include "./GranularSynth.hpp"
#include "./Util.hpp"
#include <iostream>
#include <cmath>

GranularSynth::GranularSynth(int table_size, float fs) :
win(Window::Hann,64), N{table_size}, fs{fs} {
    // std::cout << "grain 1 size = " << grains[1].size;
        gainTracker.setT(1000);
        ratioTracker.setT(1000);
    grains.resize(table_size);
}

//GranularSynth::GranularSynth() :
//win(Window::Hann, 64), N{100} {
//    grains.resize(N);
//}

void GranularSynth::newGrain(float *data, float offset, float start_ix) {
    int found_grain = 0;

    // linear search for first open grain
    for (int k = 0; k < N; k++){
        if (grains[k].size > 0) {
            continue;
        }
        //std::cout << "found open grain at " << k << std::endl;
       // std::cout << "new grain start_ix = " << start_ix << std::endl;
        grains[k].size = length;
        grains[k].ix = start_ix;
        grains[k].data = data;
        grains[k].win_ix = 0;
        grains[k].offset = offset;

        if (k > maxgrain)
            maxgrain = k;
        
        found_grain = true;
        break;
    }

    if (!found_grain) {
        miss_count++;
    }

    for (int k = maxgrain; k > 0; k--){
        if (grains[k].size >= 0)
            break;
        
        maxgrain = k;
    }
}

float GranularSynth::setGain(float gain) {
    this->gain = gain;
    return gain;
}

float GranularSynth::setRatio(float ratio) {
    this->ratio = ratio;
    return ratio;
}

float GranularSynth::setT(float t) {
    this->T = t;
    return T;
}

float GranularSynth::setPan(float pan) {
    this->pan = pan;
    return pan;
}

void GranularSynth::setGrainSource(float *data, float offset, float length) {
    this->source = data;
    this->length = length;
    this->offset = offset;
}

float* GranularSynth::getGrainSource() {
    return source;
}

float GranularSynth::setVibratoAmpl(float ampl) {
    vib_a = ampl;
    return ampl;
}

float GranularSynth::setVibratoRate(float freq) {
    vib_f = freq;
    return freq;
}

float GranularSynth::synthesizeOne() {
    
    theta += vib_f * M_2_PI/fs;
    if (theta > M_2_PI)
        theta = 0.f;
    
    if ( enable && (--nextgrain < 0)) {
        newGrain(this->source,offset,-nextgrain);
        nextgrain += T + vib_a*sinf(theta);
    }
    float smooth_gain = gainTracker.compute(gain);
    float smooth_ratio = ratioTracker.compute(ratio);
    
    float sum = 0.0;
    for (int ix = 0; ix <= maxgrain; ix++) // look for active grains
    {         
        // if this grain has been "triggered", it's size is > 0  
        if (grains[ix].size <= 0) continue;
        float u = valueAtIndexInterp(grains[ix].data, grains[ix].ix + grains[ix].offset) * smooth_gain;
        float f = grains[ix].win_ix/(grains[ix].size-1);
        float w = win_enable ? win.value(f) : 1.0;  
        u *= w;
        //std::cerr << u << std::endl;     

        sum += u;
        grains[ix].ix += smooth_ratio;
        grains[ix].win_ix += smooth_ratio;
        
        if (grains[ix].ix >= grains[ix].size){
            grains[ix].size = -1; // declare this grain open
        }
    }
    return sum;
}

void GranularSynth::synthesize(float *out, int n) {
    for (int k = 0; k < n; k++){ // for each sample
        *(out + k) = this->synthesizeOne();
    }
}
