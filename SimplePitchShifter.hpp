#ifndef _simpleshifter_
#define _simpleshifter_

#include <cstdio>
#include <string>
#include "CircularAudioBuffer.hpp"
#include "Window.hpp"

class SimplePitchShifter {
public:
    SimplePitchShifter(CircularAudioBuffer &b, Window &w, float MaxT);
    ~SimplePitchShifter();

    void setPeriod(float T);
    void setRatio(float r);
    float computeOne();
    void computeReplace(float * output, int N);
    void computeAdd(float * output, int N);
    
private:
    float ix1, ix2, T, ratio;
    float xfade_ix, xfade_n;
    float maxT;
    CircularAudioBuffer& b;
    Window& w;
};

#endif