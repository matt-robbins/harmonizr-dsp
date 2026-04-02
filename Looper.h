//
//  Looper.h
//  iOSHarmonizerFramework
//
//  Created by Matthew E Robbins on 4/18/23.
//

#ifndef Looper_h
#define Looper_h

#include <stdlib.h>
#include <vector>
#include "Window.hpp"

class Looper {
public:
    enum loopMode {
        LoopStopped=0,
        LoopRec,
        LoopPlay,
        LoopPlayRec,
        LoopPause
    };
    Looper(int channels=1, int L=44100, int xfn=100, int wlen=200);
    ~Looper();

    // run on an array of channel data buffers
    int compute(float * buf[], int offset, int count);
    float position();
    loopMode setMode(loopMode mode);

    loopMode loop_mode;

private:
    Window win;
    int ix; // current position
    int L; // total length
    int n; // current loop length
    int xfn; // cross fade duration in samples
    int xf; // current x-fade position
    int channels;
    std::vector<std::vector<float>> buffers; //one per channel
};


#endif /* Looper_h */
