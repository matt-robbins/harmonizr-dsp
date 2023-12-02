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
    ButterworthFilter() {}
    
    void compute(float * out, float * in, int n);
    float compute_one(float in);
    void clear_bad();
    
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
    int ix = 0;
};



#endif /* ButterworthFilter_h */

