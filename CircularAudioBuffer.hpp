#ifndef _circbuf_
#define _circbuf_

#include <cstdio>
#include <string>

class CircularAudioBuffer {
public:
    CircularAudioBuffer(int N);
    ~CircularAudioBuffer();

    void insertValueAtIndex(int ix, float val);
    int pushValue(float val);
    float valueAtIndexInterp(float ix);

    float wrapIndex(float ix);

    const std::string toString();

    float operator[](int idx);

    void copyRange(int ix, int N, float *out);
    int getWriteIndex();
    int getSize();

    float relIndexBehind(float ix);
    
private:
    float *data;
    float *_data;
    int N;
    int ix;
};

#endif
