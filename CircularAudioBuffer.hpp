#ifndef _circbuf_
#define _circbuf_

#include <cstdio>
#include <string>

class CircularAudioBuffer {
public:
    CircularAudioBuffer(int N);
    ~CircularAudioBuffer();

    void insertValueAtIndex(int ix, float val);
    int insertValue(float val);
    float valueAtIndexInterp(float ix);

    const std::string toString();

    float operator[](int idx);

    void copyRange(int ix, int N, float *out);
    
private:
    float *data;
    int N;
    int ix;
};

#endif