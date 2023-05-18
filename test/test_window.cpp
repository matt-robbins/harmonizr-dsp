#include "../Looper.h"
#include "common.h"

namespace {

class LooperTest : public testing::Test {
  protected:
    void SetUp() override {
      
      buf = new float*[nch];
      for (int j = 0; j < nch; j++) {
        buf[j] = new float[N]{};
      }

      int v = 1;
      for (int k = 0; k < N/2; k++){
          buf[0][k] = k - N/8;
          buf[1][k] = k - N/8;
          v *= -1;
      }
      for (int k = N/2; k < 3*N/4; k++){
          buf[0][k] = 1.0;
          buf[1][k] = 1.0;
          v *= -1;
      }
    }

    void TearDown() override {
      for (int j = 0; j < nch; j++) {
        delete(buf[j]);
      }
      delete(buf);
    }
    int N = 32;
    int nch = 2;
    Looper l = Looper(2,N/4,N/8,nullptr);
    
    float **buf;
};

TEST_F(LooperTest, TestLoop) {
  l.setMode(Looper::LoopRec);
  
  int n = N/4;

  l.compute(buf,0,n);
  l.setMode(Looper::LoopPlay);
  l.compute(buf,n,n);
  l.setMode(Looper::LoopPlayRec);
  l.compute(buf,2*n,n);
  l.setMode(Looper::LoopPause);
  l.compute(buf,3*n,n);
  for (int k = 0; k < N; k++) {
    for (int ch = 0; ch < nch; ch++) {
      std::cout << buf[ch][k] << ", ";
    }
    std::cout << "\n";
  }
}

TEST(TestWindow, TestCenter) {
  int n_instances = 0;

  // Expect our window to return a value of 1 for the center position
  Window w = Window(Window::Hann, 64);
  float d = 1-w.value(0.5);

  EXPECT_LT(d,0.00001) << ERR_PREFIX << "Fail!";
  
  d = 0.5-w.value(0.25);

  EXPECT_LT(d,0.00001) << ERR_PREFIX << "Fail!";

  EXPECT_EQ(0, w.value(0)) << ERR_PREFIX << "Fail!";
  EXPECT_EQ(0, w.value(1)) << ERR_PREFIX << "Fail!";
}

int main(int argc, char** argv) {
  testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();

}

} // namespace