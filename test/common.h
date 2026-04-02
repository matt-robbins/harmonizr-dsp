// Common for all tests
#ifndef TEST_INC_COMMON_H_
#define TEST_INC_COMMON_H_

#include "gtest/gtest.h"

const char *ERR_PREFIX = ">> ";     // For printed error descriptions
const char *NOTE_PREFIX = "  $$ ";  // For printed notes/comments

class ScopeTimer {
  std::chrono::time_point<std::chrono::steady_clock> start, end;
  std::chrono::duration<float> duration;
  std::string m_name;
  public:
  ScopeTimer(std::string&& name) {
      m_name = name;
      start = std::chrono::high_resolution_clock::now();
  }
  ~ScopeTimer() {
      end = std::chrono::high_resolution_clock::now();
      duration = end - start;
      float ms = duration.count() * 1000.0f;
      std::cout << m_name << ": " << ms << "ms" << std::endl;
  }
};
#endif  // TEST_INC_COMMON_H_