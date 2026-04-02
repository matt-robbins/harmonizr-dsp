#include "../GranularSynth.hpp"
#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"

namespace {

TEST(TestSynth, TestSynth) {
    int n_instances = 0;

    GranularSynth synth = GranularSynth(5);
    synth.win_enable = true;

    float buf[32]{};
    float data[16] = {0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0};

    synth.setGrainSource(data+3,-2,13);
    synth.T = 4.25;

    //std::memset(buf, 0, 32*sizeof(float));

    synth.synthesize(buf, 32);    

    EXPECT_LT(fabs(buf[5]-1.0),0.001);
    EXPECT_LT(fabs(buf[13]-0.5529),0.001);
    
    for (int j = 0; j < 32; j++) {
        std::cerr << buf[j] << std::endl;
    }
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}