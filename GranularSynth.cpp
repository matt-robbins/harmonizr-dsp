#include "./GranularSynth.hpp"
#include <iostream>

GranularSynth::GranularSynth(int table_size) : N{table_size} {
   // std::cout << "grain 1 size = " << grains[1].size;
   grains.resize(table_size);
}

GranularSynth::~GranularSynth() {
}

void GranularSynth::newGrain(
    CircularAudioBuffer *source, float start_ix, int size, float gain, float ratio, float pan) {
    int found_grain = 0;

    // linear search for first open grain
    for (int k = 0; k < N; k++){
        if (grains[k].size >= 0) {
            continue;
        }
        //std::cout << "found open grain at " << k << std::endl;
        grains[k].size = size;
        grains[k].start = start_ix;
        grains[k].ratio = ratio;
        grains[k].gain = gain;
        grains[k].pan = pan;

        grains[k].ix = 0;
        grains[k].b = source;

        if (k > maxgrain)
            maxgrain = k;
        
        found_grain = true;
        break;
    }

    for (int k = maxgrain; k > 0; k--){
        if (grains[k].size >= 0)
            break;
        
        maxgrain = k;
    }
}

float GranularSynth::synthesizeOne() {
    float sum = 0.0;
    for (int ix = 0; ix <= maxgrain; ix++) // look for active grains
    {            
        if (grains[ix].size > 0){ // if this grain has been "triggered", it's size is > 0   
            
            std::cout << "found active grain at index " << ix << std::endl;
            float u = grains[ix].b->valueAtIndexInterp(grains[ix].start + grains[ix].ix);
            float f = grains[ix].ix/(grains[ix].size-1);            
            u = u * win.value(f);
            
            sum += u;
            grains[ix].ix += grains[ix].ratio;
            
            if (grains[ix].ix > grains[ix].size){
                grains[ix].size = -1; // declare this grain open
            }
        }
    }
    return sum;
}

void GranularSynth::synthesize(float *out, int n) {
    for (int k = 0; k < n; k++){ // for each sample
        *(out + k) = this->synthesizeOne();
    }
}