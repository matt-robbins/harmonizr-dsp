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

    float computeOne();
    void computeReplace(float * output, int N);
    void computeAdd(float * output, int N);
    float ratio;
    float T;
private:
    float ix1, ix2;
    float xfade_ix, xfade_n;
    float maxT;
    CircularAudioBuffer& b;
    Window& w;
};

#endif
