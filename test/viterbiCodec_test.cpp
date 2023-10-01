#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <iostream>
#include "../src/feedForwardTrellis.h"
#include "../src/viterbiCodec.h"
#include "../src/minHeap.h"

namespace {
std::vector<double> addNoise(std::vector<int> modulated_signal, double SNR) {
  std::random_device rd;
  std::mt19937 noise_gen( 82 );
  std::vector<double> noisyMsg;

  /*
  // the below lines break the noise generator for some compilers
  std::random_device rd;
  std::default_random_engine generator;
  generator.seed( rd() ); //Now this is seeded differently each time.
  */
  double variance = pow(10.0, -SNR / 10.0);
  double sigma = sqrt(variance);
  std::normal_distribution<double> distribution(0.0, sigma);

  for (int i = 0; i < modulated_signal.size(); i++) {
    noisyMsg.push_back(modulated_signal[i] + distribution(noise_gen));
  }
  return noisyMsg;
}
}

TEST(ViterbiCodec, listDecoding) {
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 3;
  code.list_size = 4;
  code.generator_poly = {13, 17};

  ViterbiCodec codec(code);
  std::vector<int> noiseless_signal = {-1, -1, -1, 1, -1, 1, -1, -1};
  std::vector<double> noisy_signal = addNoise(noiseless_signal, 0.0);
  std::vector<MessageInformation> output = codec.listViterbiDecoding(noisy_signal);
  
  std::vector<int> path_1 = output[0].path;
  std::vector<int> message_1 = output[0].message;
  std::vector<int> path_2 = output[1].path;
  std::vector<int> message_2 = output[1].message;
  std::vector<int> path_3 = output[2].path;
  std::vector<int> message_3 = output[2].message;
  std::vector<int> path_4 = output[3].path;
  std::vector<int> message_4 = output[3].message;
  std::cout << "printing path 1: ";
  CodecUtils::print(path_1);
  CodecUtils::print(message_1);

  EXPECT_EQ(path_1.front(), 0);
  EXPECT_EQ(path_1.back(), 5);

  std::cout << "printing path 2: ";
  CodecUtils::print(path_2);
  EXPECT_EQ(message_2.size(), 4);
  CodecUtils::print(message_2); 
  
  EXPECT_EQ(path_2.front(), 0);
  EXPECT_EQ(path_2.back(), 1);

  std::cout << "printing path 3: ";
  CodecUtils::print(path_3);

  CodecUtils::print(message_3);

  std::cout << "printing path 4: ";
  CodecUtils::print(path_4);

  CodecUtils::print(message_4);

}

TEST(ViterbiCodec, CRCTest) {

  std::vector<int> input = {1, 1, 0, 0, 1, 0, 1, 0, 1};
  int crc_dec = 21;
  int crc_length = 5;
  std::vector<int> crc_message = CRC::calculateCRC(input, crc_dec, crc_length);
  CodecUtils::print(crc_message);

  bool all_zero = CRC::checkCRC(crc_message, crc_dec, crc_length);
  all_zero = CRC::checkCRC({1, 1, 0, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0}, crc_dec, crc_length);
  EXPECT_EQ(all_zero, false);

  

}




int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}