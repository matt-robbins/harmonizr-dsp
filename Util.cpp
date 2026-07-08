//
//  Util.cpp
//  Harmonizer
//
//  Created by Matthew E Robbins on 7/8/26.
//
#include "Util.hpp"

float midi_note_to_T(float nn, float fs) {
    float f = 440.f * powf(2.f, (nn - 69.f)/12.f);
    float T = fs/f;
    return T;
}
