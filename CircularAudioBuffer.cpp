#include "CircularAudioBuffer.hpp"
#include "Util.hpp"
#include <cstring>
#include <algorithm>
#include <iostream>

CircularAudioBuffer::CircularAudioBuffer(int log2N) : log2N{log2N} {
    if (log2N > 16) {
        std::cerr << "sizes bigger than 2^16 make no sense.\n";
        log2N = 16;
    }
    if (log2N < 4) {
        std::cerr << "sizes smaller than 2^4 make no sense.\n";
        log2N = 4;
    }
    //std::cerr << "constructing CircBuf with size " << log2N << std::endl;
    N = 0x1 << log2N;
    mask = N-1;
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
        if (i == ix || i == ix+N) {
            out += "-->";
        }
        out += std::to_string(data[i]) + "\n";
    }
    return out;
}

float CircularAudioBuffer::valueAtIndexInterp(float ix) {

    int iix = (int) floorf(ix); 
    float frac = ix - iix;
    iix = iix&mask;
    if ((iix + 3) > N) iix -= N;

    return cubic (data + iix - 1, frac); // cubic() needs one index behind
}

float CircularAudioBuffer::operator[](int ix){
    return data[ix&mask];
}

void CircularAudioBuffer::copyRange(int relix, int n, float * out) {

    if (relix < 1 || relix > N) {
        std::cerr << "relix = " << relix << "\n";
        throw std::invalid_argument("relative index must be less than N samples behind");
    }

    if (n > relix) {
        throw std::invalid_argument("input size n must be smaller than buffer size N");
    }

    std::memcpy(out, data + ix-relix, n * sizeof(float));
}

// return a raw pointer to a congiguous segment of the buffer,
// guaranteed to be N long. Offset is absolute.
float *CircularAudioBuffer::getContiguous(int offset) {
    if (offset > 0) offset -= N;
    if (offset < -N) offset += N;
    
    return data + offset;
}

// return a raw pointer to a congiguous segment of the buffer,
// guaranteed to be N long. Offset is relative.
float *CircularAudioBuffer::getContiguousRelative(int offset) {
    offset += ix;
    if (offset > 0) offset -= N;
    if (offset < -N) offset += N;
    
    return data + offset;
}

int CircularAudioBuffer::pushValue(float val){
    int r = ix;
    this->insertValueAtIndex(ix,val);
    if (++ix >= N)
        ix = 0;

    return r;
}

int CircularAudioBuffer::pushData(float *in, int n) {
    //int r = ix;
    if (n > N){
        return -1;
    }
    int rp = std::min(n, N-ix);
    if (rp) {
        std::copy(in,in+rp,data+ix);
        std::copy(in,in+rp,_data+ix);
    }

    ix += rp;

    int rem = n-rp;
    if (rem) {
        std::copy(in+rp,in+n,data);
        std::copy(in+rp,in+n,_data);
        ix = rem;
    }
    return n;
}


void CircularAudioBuffer::insertValueAtIndex(const int idx, float val) {
    data[idx] = _data[idx] = val;
    return;
}
