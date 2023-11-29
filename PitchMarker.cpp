#include "PitchMarker.hpp"
#include <iostream>

PitchMarker::PitchMarker(CircularAudioBuffer &b, float maxT) : b{b}, maxT{maxT} {
    pix = b.getWriteIndex() - maxT/2;
}
PitchMarker::~PitchMarker() { }

float PitchMarker::getMark() {
    return pix;
}

int PitchMarker::findMarkD(float T, float frac) {
    int wp = b.getWriteIndex();

    // this condition is true if and only if b has circled around to 0. 
    if (pix > wp) {
        pix -= b.getSize();
    }
    
    // return without finding new mark if it hasn't been "long enough".
    if (wp - pix < (maxT/2 + T*1.25))
        return pix;

    // swap marks
    old_pix = pix; 
    pix += T;

    // Jump forward to the latest if there is too much space
    while (wp - pix > (maxT/2 + 2*T)) {
        std::cout << "jumping..." << std::endl;
        pix += T;
    }

    int srch_rng = (int) T * frac;

    // get minimum
    float min = HUGE_VALF;
    for (float k = -srch_rng; k <= srch_rng; k++){
        float val = b[pix+k];
        if (val < min)
            min = val;
    }
    
    // get sum, subtracting off min
    float sum = 0;
    for (float k = -srch_rng; k <= srch_rng; k++){
        float add = b[pix+k] - min;
        sum += add;
    }

    // find median, as in a distribution
    float csum = 0, last_csum = 0;
    float median = 0;
    for (float k = -srch_rng; k <= srch_rng; k++){
        csum += (b[pix+k] - min);
        if (csum <= sum/2){
            last_csum = csum;
            continue;
        }

        // linearly interpolate
        float frac = ((sum/2)-last_csum)/(csum-last_csum);
        median = k-1+frac; 
        break; 
    }

    std::cout << "median: " << median << std::endl;
    if (sum > 0)
        pix += median;
        
    return pix;
}

float PitchMarker::findMark(float T, float res) {

    int wp = b.getWriteIndex();

    // this condition is true if and only if b has circled around to 0. 
    if (pix > wp) {
        pix -= b.getSize();
    }
    
    // return without finding new mark if it hasn't been "long enough".
    if (wp - pix < (maxT/2 + T*1.25))
        return pix;

    // swap marks
    old_pix = pix; 
    pix += T;

    // Jump forward to the latest if there is too much space
    while (wp - pix > (maxT/2 + 2*T)) {
        std::cout << "jumping..." << std::endl;
        pix += T;
    }
        

    // find center of mass around next pitch mark
    float srch_rng = T/2.0;
    float srch_inc = srch_rng/res;

    // get minimum
    float min = HUGE_VALF;
    for (float k = -srch_rng; k <= srch_rng; k+=srch_inc){
        float val = b.valueAtIndexInterp(pix+k);
        if (val < min)
            min = val;
    }
    
    // get sum, subtracting off min
    float sum = 0;
    for (float k = -srch_rng; k <= srch_rng; k+=srch_inc){
        float add = (b.valueAtIndexInterp(pix+k) - min);
        sum += add;
    }

    // find median, as in a distribution
    float csum = 0, last_csum = 0;
    float median = 0;
    for (float k = -srch_rng; k <= srch_rng; k+=srch_inc){
        csum += (b.valueAtIndexInterp(pix+k) - min);
        if (csum <= sum/2){
            last_csum = csum;
            continue;
        }

        // linearly interpolate
        float frac = ((sum/2)-last_csum)/(csum-last_csum);
        median = k-1+frac; 
        break; 
    }

    std::cout << "median: " << median * srch_inc << std::endl;
    
    if (sum > 0)
        pix += median * srch_inc;

    return pix;
}