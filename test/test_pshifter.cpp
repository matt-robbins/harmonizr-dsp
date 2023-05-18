#include "../SimplePitchShifter.hpp"
#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"
#include <cmath>

namespace {

TEST(TestCbuf, TestCenter) {
    int n_instances = 0;

    int N = 32;
    int T = 12;
    Window w = Window(Window::Hann, N);
    CircularAudioBuffer b = CircularAudioBuffer(N);

    SimplePitchShifter p = SimplePitchShifter(b,w,16);
    p.setPeriod(T);
    p.setRatio(12.0/23.0f);

    for (int k = 0; k < N; k++) {
        b.insertValue(sinf(M_PI*2*k/T));
    }

    for (int k = 0; k < N; k++) {
        std::cout << p.computeOne() << std::endl;
        //b.insertValue(sinf(M_PI*2*(k+N)/T));
    }

}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}