#include "../PitchEstimator.hpp"
#include "../PitchMarker.hpp"
#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"
#include <cmath>

namespace {


TEST(TestPitchEst, TestYIN) {
    int n_instances = 0;

    int l2n = 11;
    int N = 0x1 << l2n;
    CircularAudioBuffer b = CircularAudioBuffer(l2n);

    std::cout << "buflen = " << b.getSize() << "\n";

    int maxT = 400;

    PitchEstimatorYIN py = PitchEstimatorYIN(maxT,l2n,0.3,1);
    PitchEstimator &p = py;

    int count = 0;
    for (float trueT = 33.1; trueT < maxT; trueT += 7.23476) {
        int offset = 0;
        for (int k = 0; k < N + N/8 + 7; k++) {
            float val = cos(k * 2*M_PI / trueT) + sin(k * 4*M_PI/trueT);
            b.pushValue(val);
        }

        float estT = p.estimate(b);
        for (int k = 0; k < 1000; k++){
            estT = p.estimate(b);
            count++;
        }

        EXPECT_LT(fabs(estT - trueT),1.0);
       // std::cout << trueT << ": " << estT << ":" << fabs(estT - trueT) << std::endl;  
    }    
    std::cout << count << "x" << std::endl;
}

TEST(TestPitchEst, TestMark) {
    int l2n = 10;
    int N = 0x1 << l2n;
    CircularAudioBuffer b = CircularAudioBuffer(l2n);

    int maxT = 400;

    PitchMarker m = PitchMarker(b,(float)maxT);

    float trueT = 200;
    float phase = 0.0;
    float sign = 1;
    float val = 0;

    float data[5] = {5,5,10,6,5};

    float median = median_idx(data,5);

    //std::cout << "pitch mark: " << m.getMark() << std::endl;

    float old_pix = 0;
    for (int k = 0; k < 100; k++) {
        for (int k = 0; k <= 128; k++) {
            float val = cos(phase * 2*M_PI / trueT);
            // triangle wave
            // val += sign * 4/trueT;
            // if (std::abs(val) >= 1) {
            //     sign = -sign;
            // }
            //std::cerr << "val=" << val << "\n";
            b.pushValue(val);
            phase += 1;
        }

        // give slightly wrong estimate for T
        float pix = m.findMark(trueT-10, 0.5);

        float diff = pix - old_pix;
        if (diff < 0)
            diff += N;

        std::cerr << " ***** pix moved forward by: " << diff << "\n";

       // EXPECT_LT(fabs(pix-old_pix - trueT),5);
        old_pix = pix;

       // std::cout << b.getWriteIndex() - pix << " samples back, " << pix << std::endl;
    }
}


int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}