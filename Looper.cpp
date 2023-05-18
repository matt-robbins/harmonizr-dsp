#include "Looper.h"
#include "Window.hpp"
#include <iostream>

Looper::Looper(int channels, int L_samples, int xfade_samples, Window *w=nullptr) : 
    channels{channels}, L{L_samples}, xfn{xfade_samples}, win{w} {

    buffers = new float*[channels];
    for (int k = 0; k < channels; k++){
        buffers[k] = new float[L];
    }

    ix = 0;
    n = 0;
    xf = 0;
    loop_mode = LoopStopped;

    if (win == NULL) {
        dwin = win = new Window(Window::Hann, 64);
    }

    std::cout << "created Looper with " << channels << " channels, " << L << " samples long" << "\n";
}

Looper::~Looper() {
    for (int k = 0; k < channels; k++) {
        delete(buffers[k]);
    }
    delete(buffers);
    if (dwin != nullptr)
        delete(dwin);
}

int Looper::compute(float * iobuf[], int offset, int count){
    for (int chan = channels - 1; chan >= 0; chan--){
        int cix = ix;
        int cxf = xf;
        for (int frameIndex = offset; frameIndex < count+offset; ++frameIndex)
        {
            float stored = buffers[chan][cix];
            float input = iobuf[chan][frameIndex];
            float w = 0;
            if (cxf > 0){
                cxf--;
                w = win->value((float) cxf / (2*xfn));
            }
            
            switch (loop_mode)
            {
                case LoopRec:
                    if (cxf > 0) {
                        input *= (1 - w);
                    }
                    buffers[chan][cix] = input; cix++;
                    break;
                case LoopPlayRec:
                    if (cxf > 0) {
                        input *= (1 - w);
                    }
                    buffers[chan][cix] += input;
                    iobuf[chan][frameIndex] += stored; cix++;
                    break;
                case LoopPlay:
                    iobuf[chan][frameIndex] += stored;
                    if (cxf > 0) {
                        buffers[chan][cix] += input * w;
                    }
                    cix++;
                    break;
                case LoopPause:
                    if (cxf > 0) {
                        iobuf[chan][frameIndex] += stored * w;
                    }
                    break;
                default:
                    break;
            }
            
            if ((cix >= L) || ((cix >= n) && (n > 0)))
            {
                cix = 0;
            }
        }
        if (chan == 0) {
            ix = cix; xf = cxf;
        }
    }

    std::cout << "contents up to size = " << L << "\n";
    for (int k = 0; k < L; k++) {
        for (int ch = 0; ch < channels; ch++) {
            std::cout << buffers[ch][k] << ", ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    return count;
}

Looper::loopMode Looper::setMode(loopMode mode) {
    loopMode old_mode = loop_mode;
    loop_mode = mode;

    std::cout << "switching from " << old_mode << " to " << loop_mode << "\n";

    if ((old_mode == LoopRec) && (n == 0)){
        n = ix;
        xf = xfn;
    }
    if (loop_mode == LoopStopped){
        n = 0;
        ix = 0;
    }
    if (loop_mode == LoopRec || loop_mode == LoopPlayRec || loop_mode == LoopPause) {
        xf = xfn;
    }
    
    return loop_mode;
}
    
float Looper::position() {
    if (n <= 0){
        return (float) ix / (float) L;
    }
    return (float) ix / (float) n;
}