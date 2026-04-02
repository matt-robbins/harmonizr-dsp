#ifndef _PitchMarker_h
#define _PitchMarker_h
#include "CircularAudioBuffer.hpp"

class PitchMarker {
    public:
        PitchMarker(CircularAudioBuffer &b, float maxT=800);
        ~PitchMarker();
        bool findMark(float T, float frac);
        float getMark();
        float mark;
    private:
        CircularAudioBuffer &b;
        float maxT;
        float old_mark;
        float res;
};

#endif
