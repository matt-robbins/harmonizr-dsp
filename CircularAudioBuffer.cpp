#include "CircularAudioBuffer.hpp"
#include "Util.hpp"
#include <cstring>
#include <algorithm>
#include <iostream>

CircularAudioBuffer::CircularAudioBuffer(int length) : N{length} {
    data = new float[N+3]{};
    ix = 0;
}

CircularAudioBuffer::~CircularAudioBuffer() {
    delete(data);
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
float CircularAudioBuffer::relIndexAhead(float c) {
    float d = c - ix;
    while (d < 0)
        d += N;

    return d;
}

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

    if (n > N) {
        throw std::invalid_argument("input size n must be smaller than buffer size N");
    }

    int startix = relix+ix;
    while (startix < 0) startix += N;
    while (startix >= N) startix -= N;
    int nback = std::min(n,N - startix - 1);
    
    std::memcpy(out, data + startix, nback * sizeof(float));
    int rem = n - nback;
    if (rem > 0) {
        std::memcpy(out + nback, data, rem * sizeof(float));
    }
}

int CircularAudioBuffer::insertValue(float val){
    int r = ix;
    this->insertValueAtIndex(ix,val);
    if (++ix >= N)
        ix = 0;

    return r;
}

void CircularAudioBuffer::insertValueAtIndex(const int idx, float val) {
    data[idx] = val;
    if (idx < 3){
        data[N+idx] = val;
    }
    return;
}