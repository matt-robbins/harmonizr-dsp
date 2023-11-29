#ifndef _PitchMarker_h
#define _PitchMarker_h
#include "CircularAudioBuffer.hpp"

class PitchMarker {
    public:
        PitchMarker(CircularAudioBuffer &b, float maxT=800);
        ~PitchMarker();
        float findMark(float T, float res);
        int findMarkD(float T, float frac);
        float getMark();
    private:
        CircularAudioBuffer &b;
        float maxT;
        float pix, old_pix;
        float res;
};

#endif