#include "CircularAudioBuffer.hpp"
#include "Util.hpp"
#include <cstring>
#include <algorithm>
#include <iostream>

CircularAudioBuffer::CircularAudioBuffer(int length) : N{length} {
    _data = new float[2*N]{}; //alias data back one full length
    data = _data + N;
    ix = 0;
}

CircularAudioBuffer::~CircularAudioBuffer() {
    delete(_data);
}

int CircularAudioBuffer::getWriteIndex() {
    return ix;
}

int CircularAudioBuffer::getSize() {
    return N;
}

float CircularAudioBuffer::relIndexBehind(float c) {
    float d = ix - c;
    while (d < 0)
        d += N;

    return d;
}
// float CircularAudioBuffer::relIndexAhead(float c) {
//     float d = c - ix;
//     while (d < 0)
//         d += N;

//     return d;
// }

float CircularAudioBuffer::wrapIndex(float cx) {
    if (cx < 0) cx += N;
    if (cx >= N) cx -= N;
    return cx;
}

const std::string CircularAudioBuffer::toString() {
    std::string out = "";

    for (int i = 0; i < N; i++) {
        if (i == ix) {
            out += "-->";
        }
        out += std::to_string(data[i]) + "\n";
    }
    return out;
}

float CircularAudioBuffer::valueAtIndexInterp(float ix) {
    float fi = ix - 1; // cubic() needs one index behind

    while (fi >= N)
        fi -= N;
    while (fi < 0)
        fi += N;
    
    int i = (int) fi;
    float u = cubic (data + i, fi - i);
    return u;
}

float CircularAudioBuffer::operator[](int ix){
    return data[ix];
}

void CircularAudioBuffer::copyRange(int relix, int n, float * out) {

    if (relix < 1 || relix > N) {
        throw std::invalid_argument("relative index must be less than N samples behind");
    }

    if (n > relix) {
        throw std::invalid_argument("input size n must be smaller than buffer size N");
    }

    std::memcpy(out, data + ix-relix, n * sizeof(float));
}

int CircularAudioBuffer::pushValue(float val){
    int r = ix;
    this->insertValueAtIndex(ix,val);
    if (++ix >= N)
        ix = 0;

    return r;
}

void CircularAudioBuffer::insertValueAtIndex(const int idx, float val) {
    data[idx] = _data[idx] = val;
    return;
}
