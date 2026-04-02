#ifndef _gransynth_
#define _gransynth_
#include "CircularAudioBuffer.hpp"
#include "Window.hpp"
#include <vector>

struct Grain
{
    float size = -1.0;
    float ix = 0.0;
    float win_ix = 0;
    float * data = nullptr;
    float offset = 0;
};

class GranularSynth {
public:
    GranularSynth(int table_size);
    //~GranularSynth();

    void newGrain(float * data, float offset, float start_ix);
    float synthesizeOne();
    void synthesize(float *out, int n);
    void setGrainSource(float *data, float offset, float length);
    float* getGrainSource();
    float setVibratoRate(float freq);
    float setVibratoAmpl(float amp);
    bool enable = true;
    bool win_enable = true;
    
    float T = 0;
    float gain = 1.0;
    float ratio = 1.0;
    
private:

    int N;
    std::vector<Grain> grains;
    int maxgrain = 0;
    float nextgrain = 0.0;
    float * source = nullptr;
    float length = 0.0;
    float offset = 0.0;
    float fs;
    
    float vib_f = 0.f;
    float vib_a = 0.f;
    float theta = 0.f;

    Window win;
};

#endif
