//
//  Looper.h
//  iOSHarmonizerFramework
//
//  Created by Matthew E Robbins on 4/18/23.
//

#ifndef Looper_h
#define Looper_h

#include <stdlib.h>
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
    Looper(int channels, int L, int xfn, Window *w);
    ~Looper();

    // run on an array of channel data buffers
    int compute(float * buf[], int offset, int count);
    float position();
    loopMode setMode(loopMode mode);

    loopMode loop_mode;

private:
    Window *win;
    Window *dwin;
    int ix; // current position
    int L; // total length
    int n; // current loop length
    int xfn; // cross fade duration in samples
    int xf; // current x-fade position
    int channels;
    float ** buffers; //one per channel
};


#endif /* Looper_h */
