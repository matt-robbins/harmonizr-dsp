#include "PitchMarker.hpp"
#include <iostream>
#include <algorithm>
#include <numeric>
#include "Util.hpp"

PitchMarker::PitchMarker(CircularAudioBuffer &b, float maxT) : b{b}, maxT{maxT} {
    //std::cerr << "buffer size " << b.getSize() << std::endl;
    //std::cerr << "maxT" << maxT << std::endl;
    mark = b.getWriteIndex() - maxT/2;
    //std::cerr << "mark " << mark << std::endl;
}
PitchMarker::~PitchMarker() { }

float PitchMarker::getMark() {
    return mark;
}

bool PitchMarker::findMark(float T, float frac) {
    int wp = b.getWriteIndex();

    // this condition is true if and only if b has circled around to 0. 
    if (mark > wp) {
        //std::cerr << "jumping back around\n";
        mark -= b.getSize();
    }

    // return without finding new mark if it hasn't been "long enough".
    if (wp - mark < (maxT + T*1.25)){
        //std::cout << "wp-mark <" << (maxT + T*1.25) << std::endl;
        //std::cerr << "not finding! wp = " << wp << " pix = " << pix << "\n";
        return false;
    }
    
    // swap marks
    old_mark = mark;
    mark += T;

    // Jump forward to the latest if there is too much space
    while (wp - mark > (3*maxT/2)) {
        //std::cerr << "jumping..." << std::endl;
        mark += T;
    }

    int srch_rng = (int) floorf(T * frac);
    

    int vlen = srch_rng * 2 + 1;
    float * data = b.getContiguous(floorf(mark - srch_rng));
    float shunt_frac = mark - floorf(mark);

    // find median, as in a distribution
    float median = median_idx(data,vlen) - srch_rng - shunt_frac;
   // std::cout << "median=" << median << "\n";
    mark += median;
    return true;
}
