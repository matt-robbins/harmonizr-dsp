#include "PitchMarker.hpp"

PitchMarker::PitchMarker(CircularAudioBuffer &b, float maxT) : b{b}, maxT{maxT} {
    pix = b.getWriteIndex() - maxT/2;
}
PitchMarker::~PitchMarker() { }

float PitchMarker::getMark() {
    return pix;
}

float PitchMarker::findMark(float T, float res) {

    int wp = b.getWriteIndex();

    // this condition is true if and only if b has circled around to 0. 
    if (pix < wp) {
        pix -= b.getSize();
    }
    
    // return without finding new mark if it hasn't been "long enough".
    if (wp - pix < (maxT/2 + T*1.25))
        return pix;

    // swap marks
    old_pix = pix;  
    pix += T;

    // find center of mass around next pitch mark
    float srch_rng = T/2.0;
    float srch_inc = srch_rng/res;

    // get minimum
    float min = HUGE_VALF;
    for (float k = -srch_rng; k <= srch_rng; k+=srch_inc){
        if (float val = b.valueAtIndexInterp(pix+k) < min)
            min = val;
    }
    
    // get sum, subtracting off min
    float sum = 0;
    for (float k = -srch_rng; k <= srch_rng; k+=srch_inc){
        sum += (b.valueAtIndexInterp(pix+k) - min);
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
    
    if (sum == 0)
        return;

    pix += median;
}