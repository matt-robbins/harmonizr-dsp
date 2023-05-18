#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"

namespace {

TEST(TestCbuf, TestCenter) {
    int n_instances = 0;

    CircularAudioBuffer b = CircularAudioBuffer(16);
    b.insertValue(1.0);
    b.insertValue(1.0);
    b.insertValue(0);
    b.insertValue(2.0); b.insertValue(3.0);
    // [1, 1, 0, 2, 3, *0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    float res = b.valueAtIndexInterp(-0.5);

    float buf[4] = {0, 0, 1, 1};
    float r2 = cubic(buf, 0.5);

    float buf2[16];
    b.copyRange(-5,8,buf2);

    std::cout << b.toString();

    EXPECT_EQ(b[-5],0);

    EXPECT_EQ(buf2[3],2.0);
    EXPECT_EQ(buf2[4],3.0);

    EXPECT_EQ(res, r2);
}

TEST(TestCbuf, TestBoost) {
    
}

int main(int argc, char** argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

}