#include "../Harmonizer.hpp"
#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"

namespace {

TEST(TestHarmonizer, TestHarmonizer) {
    int n_instances = 0;

    Harmonizer h = Harmonizer(7,1,30);
    h.setPitchEstPeriod(8);

    int per = 15;
    float *data = (float *)calloc(per,sizeof(float));
    float *out = (float *)calloc(per,sizeof(float));

    for (int k = 0; k < per; k++) {
        data[k] = sinf((float)k/2.0)/(k+1);
    }
    
    h.setVoiceT(0,20);

    for (int k = 0; k < 30; k++){
        h.compute(&data,&out,1,per);
        if (k < 5)
            continue;

        for (int k = 0; k < per; k++){
            std::cerr << out[k] << "\n";
        }
    }

    free(data);
    free(out);
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}