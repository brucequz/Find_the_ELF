#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "../src/feedForwardTrellis.h"
#include "../src/viterbiCodec.h"
#include "../src/minHeap.h"

TEST(ViterbiCodec, addNoiseTest) {
  std::vector<int> test_1 = {0, 1, 1, 0, 1, 0};
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);

}






int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}