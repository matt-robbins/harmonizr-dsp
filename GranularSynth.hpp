#ifndef _gransynth_
#define _gransynth_
#include "CircularAudioBuffer.hpp"
#include "Window.hpp"

typedef struct Grain
{
    float size, start, ix, ratio, gain, pan;

    CircularAudioBuffer *b = nullptr;
    Grain(){
        size = -1.0;
        start = -1.0;
        ix = 0.0;
        ratio = 1.0;
        gain = 1.0;
        pan = 0.0;
        b = nullptr;
    };
} Grain;

class GranularSynth {
public:
    GranularSynth(int table_size);
    ~GranularSynth();

    void newGrain(CircularAudioBuffer *b, float start_ix, int size, float gain, float ratio, float pan);
    float synthesizeOne();
    void synthesize(float *out, int n);
private:

    int N;
    Grain *grains;
    int maxgrain = 0;

    Window win = Window(Window::Hann, 64);
};

#endif