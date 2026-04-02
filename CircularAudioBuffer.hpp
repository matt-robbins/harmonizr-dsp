#ifndef _circbuf_
#define _circbuf_

#include <cstdio>
#include <string>

class CircularAudioBuffer {
public:
    CircularAudioBuffer(int log2N=11);
    ~CircularAudioBuffer();

    void insertValueAtIndex(int ix, float val);
    int pushValue(float val);
    int pushData(float *data, int N);
    float valueAtIndexInterp(float ix);
    float wrapIndex(float ix);
    const std::string toString();
    float operator[](int idx);
    void copyRange(int ix, int N, float *out);
    int getWriteIndex();
    int getSize();
    float *getContiguous(int offset);
    float *getContiguousRelative(int offset);
    float relIndexBehind(float ix);
    
private:
    float *data;
    float *_data;
    int N;
    int log2N;
    unsigned mask;
    int ix;
};

#endif
