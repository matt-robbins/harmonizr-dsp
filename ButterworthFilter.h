//
//  ButterworthFilter.h
//  iOSHarmonizerFramework
//
//  Created by Matthew E Robbins on 4/13/23.
//

#ifndef ButterworthFilter_h
#define ButterworthFilter_h

/*
    ButterworthFilter
    keeps track of state and computes via direct form I
 */

class ButterworthFilter {
public:
    ButterworthFilter() {
//        N = sizeof(b) / sizeof(float);
//        x = (float *) calloc(N, sizeof(float));
//        y = (float *) calloc(N, sizeof(float));
        ix = 0;
    }
    void compute(float * out, float * in, int n)
    {
        for (int k = 0; k < n; k++)
        {
            out[k] = compute_one(in[k]);
        }
    }
    
    inline float compute_one(float in)
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
    
    void clear_bad()
    {
        for (int k = 0; k < N; k++)
        {
            float absx = fabs(x[k]);

            if (absx < 1e-15 || absx > 1e15) {
                x[k] = 0;
            }
            
            absx = fabs(y[k]);
            
            if (absx < 1e-15 || absx > 1e15) {
                y[k] = 0;
            }
        }
    }
    
private:
    int N = 5;

    float a[5] = {1.000000000000000e+00,
        -3.847574620836099e+00,
         5.555010102851863e+00,
        -3.567181135121114e+00,
         8.597462655075311e-01};

    float b[5] = {2.651890014535641e-03,
         0.000000000000000e+00,
        -5.303780029071281e-03,
         0.000000000000000e+00,
         2.651890014535641e-03};
    float x[5] = {0, 0, 0, 0, 0};
    float y[5] = {0, 0, 0, 0, 0};
    int ix;
};



#endif /* ButterworthFilter_h */

