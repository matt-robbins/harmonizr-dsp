#include "SimplePitchShifter.hpp"
#include <iostream>


// maxT is the target latency
SimplePitchShifter::SimplePitchShifter(CircularAudioBuffer& b, Window& w, float maxT) : 
    b{b}, w{w}, maxT{maxT} 
{
    T = 400;
    ratio = 1.0;
    ix1 = 0.0f;
    ix2 = 0.0f;
    xfade_ix = 0;
    xfade_n = 0;
}

SimplePitchShifter::~SimplePitchShifter() {}

void SimplePitchShifter::setPeriod(float t) {
    T = t;
}

void SimplePitchShifter::setRatio(float r) {
    ratio = r;
}

float SimplePitchShifter::computeOne() {
    ix1 += ratio; ix2 += ratio;
    float d = b.relIndexBehind(ix1);

    // if the index gets one period behind or in front of our target latency (maxT)
    // move it by one period and crossfade

    if (d > maxT + T) {
        ix2 = ix1; ix1 += T;
        xfade_ix = xfade_n = (int) T/4;
    }
    if (d < maxT - T) {
        ix2 = ix1; ix1 -= T;
        xfade_ix = xfade_n = (int) T/4;
    }

    ix1 = b.wrapIndex(ix1);
    ix2 = b.wrapIndex(ix2);

    float val;
    if (xfade_ix > 0) {
        float mix = w.value(0.5 * (xfade_ix - 0.5) / xfade_n);
        val = mix * b.valueAtIndexInterp(ix2) + (1-mix) * b.valueAtIndexInterp(ix1);
        xfade_ix--;
    }
    else {
        val = b.valueAtIndexInterp(ix1);
    }
    
    return val;
}

void SimplePitchShifter::computeAdd(float * out, int N) {
    float d;

    for (int k = 0; k < N; k++) {
        out[k] += computeOne();
    }
}

void SimplePitchShifter::computeReplace(float * out, int N) {
    float d;

    for (int k = 0; k < N; k++) {
        out[k] = computeOne();
    }
}