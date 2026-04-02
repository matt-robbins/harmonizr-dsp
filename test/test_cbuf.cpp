#include "../CircularAudioBuffer.hpp"
#include "../Util.hpp"
#include "common.h"

namespace {

TEST(TestCbuf, TestCenter) {
    int n_instances = 0;

    int l2size = 4;
    int size = 16;
    CircularAudioBuffer b = CircularAudioBuffer(l2size);
    b.pushValue(1.0);
    b.pushValue(1.0);
    b.pushValue(0);
    b.pushValue(2.0); b.pushValue(3.0);
    // [1, 1, 0, 2, 3, *0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

    CircularAudioBuffer b2 = CircularAudioBuffer(l2size);
    float data[7] = {1,2,3,4,5,6,7};
    b2.pushData(data,7);
    b2.pushData(data,7);
    b2.pushData(data,7);

    // [3, 4, 5, 6, 7, *6, 7, 1, 2, 3, 4, 5, 6, 7, 1, 2]

    std::cout << b2.toString();


    EXPECT_EQ(b.getSize(),size);

    float res = b.valueAtIndexInterp(-1.25);

    float buf[4] = {0, 0, 1, 1};
    float r2 = cubic(buf, 0.75);

    float buf2[16];
    b.copyRange(5,5,buf2);

    std::cout << b.toString();
    EXPECT_EQ(b.getWriteIndex(),5);

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