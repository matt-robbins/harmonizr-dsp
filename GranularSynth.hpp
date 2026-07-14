#ifndef _gransynth_
#define _gransynth_
#include "CircularAudioBuffer.hpp"
#include "Window.hpp"
#include "Util.hpp"
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
    GranularSynth(int table_size=100, float fs=44100.f);

    void newGrain(float * data, float offset, float start_ix);
    float synthesizeOne();
    void synthesize(float *out, int n);
    void setGrainSource(float *data, float offset, float length);
    float* getGrainSource();
    float setVibratoRate(float freq);
    float setVibratoAmpl(float amp);
    float setGain(float gain);
    float setRatio(float ratio);
    float setPan(float pan);
    float setT(float t);
    float setFs(float fs) {
        this->fs = fs;
        return this->fs;
    }
    
    bool enable = true;
    bool win_enable = true;
    
    float T = 0;

    
private:
    int N;
    std::vector<Grain> grains;
    IirTracker gainTracker;
    IirTracker ratioTracker;
    Window win;
    
    int maxgrain = 0;
    float nextgrain = 0.0;
    float * source = nullptr;
    float length = 0.0;
    float offset = 0.0;
    
    float vib_f = 0.f;
    float vib_a = 0.f;
    float theta = 0.f;
    
    float gain = 1.0;
    float ratio = 1.0;
    float pan = 0.0;

    unsigned int miss_count = 0;
    
protected:
    float fs = 44100.f;
    
};

#endif
