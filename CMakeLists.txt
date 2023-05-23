cmake_minimum_required(VERSION 3.20)

project(HarmonizrDsp)

include(FetchContent)
FetchContent_Declare(
  googletest
  # Specify the commit you depend on and update it regularly.
  URL https://github.com/google/googletest/archive/5376968f6948923e2411081fd9372e71a59d8e77.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# Now simply link against gtest or gtest_main as needed. Eg

set(CMAKE_CXX_STANDARD 14)
include_directories(/usr/local/include)
set(SOURCES Window.cpp Looper.cpp CircularAudioBuffer.cpp 
  GranularSynth.cpp SimplePitchShifter.cpp PitchEstimator.cpp)
add_library(harm ${SOURCES})
add_executable(test_window test/test_window.cpp)
add_executable(test_cbuf test/test_cbuf.cpp)
add_executable(test_gsynth test/test_gsynth.cpp)
add_executable(test_pshift test/test_pshifter.cpp)
add_executable(test_pitchest test/test_pitchest.cpp)
target_link_libraries(test_window gtest gtest_main harm)
target_link_libraries(test_cbuf gtest gtest_main harm)
target_link_libraries(test_gsynth gtest gtest_main harm)
target_link_libraries(test_pshift gtest gtest_main harm)
target_link_libraries(test_pitchest gtest gtest_main harm "-framework Accelerate")

add_test(NAME window COMMAND test_window)
