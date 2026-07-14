#ifndef _PitchMarker_h
#define _PitchMarker_h
#include "CircularAudioBuffer.hpp"

class PitchMarker {
public:
    PitchMarker(CircularAudioBuffer &b, float maxT=800, float search_frac=0.25);
    ~PitchMarker();
    bool findMark(float T, bool voiced);
    float getMark();
    float mark;
private:
    CircularAudioBuffer &b;
    float maxT;
    float old_mark;
    float res;
    float search_frac;
};

#endif
