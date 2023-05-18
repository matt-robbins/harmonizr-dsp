#include "../GranularSynth.hpp"
#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"

namespace {

TEST(TestCbuf, TestCenter) {
    int n_instances = 0;

    // Expect our window to return a value of 1 for the center position
    GranularSynth synth = GranularSynth(3);
    CircularAudioBuffer b = CircularAudioBuffer(32);

    int offset = 0;
    for (int k = 0; k < 32; k++) {
        float val = k/4 & 0x1 ? 1.f : -1.f;
        b.insertValueAtIndex(k,val);
    }

    float buf[16]{};

    for (int j = 0; j < 1; j++){
        std::memset(buf, 0, 16*sizeof(float));
        synth.newGrain(&b, 4, 8.0f, 1.0f, 1.0f, 0.0f); 
        synth.synthesize(buf, 4);   
        synth.newGrain(&b, 28, 8.0f, 1.0f, 1.0f, 0.0f); 
        synth.synthesize(buf+4, 12);   
    }

    for (int j = 0; j < 16; j++) {
        std::cout << buf[j] << std::endl;
    }
    std::cout << std::endl;  
    
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}