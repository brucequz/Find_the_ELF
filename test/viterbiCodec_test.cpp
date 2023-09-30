#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include "../src/feedForwardTrellis.h"
#include "../src/viterbiCodec.h"
#include "../src/minHeap.h"

TEST(ViterbiCodec, listDecoding) {
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 3;
  code.list_size = 2;
  code.generator_poly = {13, 17};

  ViterbiCodec codec(code);
  std::vector<double> noiseless_signal = {-1, -1, -1, 1, -1, 1, -1, -1};
  std::vector<MessageInformation> output = codec.listViterbiDecoding(noiseless_signal);
  
  std::vector<int> path_1 = output[0].path;
  std::vector<int> path_2 = output[1].path;
  std::cout << "printing path 1: " << std::endl;
  CodecUtils::print(path_1);

  EXPECT_EQ(path_1.front(), 0);
  EXPECT_EQ(path_1.back(), 5);

  std::cout << "printing path 2: " << std::endl;
  CodecUtils::print(path_2);
  
  EXPECT_EQ(path_2.front(), 0);
  EXPECT_EQ(path_2.back(), 1);

}






int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}