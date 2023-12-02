
#include "ButterworthFilter.hpp"
#import <vector>


void ButterworthFilter::compute(float * out, float * in, int n)
{
    for (int k = 0; k < n; k++)
    {
        out[k] = compute_one(in[k]);
    }
}

inline float ButterworthFilter::compute_one(float in)
{
    float ret = 0;
    x[ix] = in;
    y[ix] = 0;
    for (int j = 0; j < N; j++)
    {
        int xix = (ix-j+N)%N;
        y[ix] += x[xix]*b[j];
    }
    
    for (int j = 1; j < N; j++)
    {
        int yix = (ix-j+N)%N;
        y[ix] -= y[yix]*a[j];
    }
    
    ret = y[ix];
    //fprintf(stderr,"%f: %f\n", in, ret);
    
    ix++;
    if (ix >= N)
    {
        ix = 0;
    }
    return ret;
}

void ButterworthFilter::clear_bad()
{
    for (int k = 0; k < N; k++)
    {
        float absx = std::fabs(x[k]);

        if (absx < 1e-15 || absx > 1e15) {
            x[k] = 0;
        }
        
        absx = std::fabs(y[k]);
        
        if (absx < 1e-15 || absx > 1e15) {
            y[k] = 0;
        }
    }
}
