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
    CircularAudioBuffer b = CircularAudioBuffer(N);

    int maxT = 400;

    PitchEstimatorYIN py = PitchEstimatorYIN(maxT,l2n,0.3,1);
    PitchEstimator &p = py;

    int count = 0;
    for (float trueT = 33.1; trueT < maxT; trueT += 7.23476) {
        int offset = 0;
        for (int k = 0; k < N + N/8 + 7; k++) {
            float val = cos(k * 2*M_PI / trueT) + sin(k * 4*M_PI/trueT);
            b.insertValue(val);
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
    int l2n = 11;
    int N = 0x1 << l2n;
    CircularAudioBuffer b = CircularAudioBuffer(N);

    int maxT = 400;

    PitchMarker m = PitchMarker(b,(float)maxT);

    float trueT = 200;
    float phase = 0.0;

    std::cout << "pitch mark: " << m.getMark() << std::endl;

    for (int k = 0; k < 100; k++) {
        for (int k = 0; k <= trueT+30; k++) {
            float val = cos(phase * 2*M_PI / trueT) + sin(phase * 4*M_PI/trueT);
            b.insertValue(val);
            phase += 1;
        }

        // give slightly wrong estimate for T
        float pix = m.findMarkD(200, 0.25);

        std::cout << b.getWriteIndex() - pix << " samples back, " << pix << std::endl;
    }
}


int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}