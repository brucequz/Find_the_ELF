#include <gtest/gtest.h>
#include <vector>
#include <random>
#include <iostream>
#include "../include/feedForwardTrellis.h"
#include "../include/viterbiCodec.h"
#include "../include/minHeap.h"

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

std::vector<int> calculateCRC(const std::vector<int>& input, int crc_dec_, int crc_length_) {
  // generating (crc_length - 1) number of redundancy bits (crc bits)
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);

  std::vector<int> output = input;
  output.resize(input.size() + crc_length_ - 1, 0);

  // long division
  for (int i = 0; i < input.size(); ++i) {
    if (output[i] == 1) {
      std::transform(output.begin() + i, output.begin() + i + crc_length_,
                     crc_bin.begin(), output.begin() + i, CRC::binSum);
    }
  }

  std::copy(input.begin(), input.end(), output.begin());

  return output;
}

bool checkCRC(std::vector<int> demodulated, int crc_dec_, int crc_length_) {
  // check crc by dividing the demodulated signal with crc poly
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);

  for (int ii = 0; ii <= (int)demodulated.size() - crc_length_; ii++) {
    if (demodulated[ii] == 1) {
      // Note: transform doesn't include .end
      std::transform(demodulated.begin() + ii,
                     demodulated.begin() + (ii + crc_length_), crc_bin.begin(),
                     demodulated.begin() + ii, CRC::binSum);
    }
  }
  bool all_zero = std::all_of(demodulated.begin(), demodulated.end(),
                              [](int i) { return i == 0; });
  return all_zero;
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
  
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 14;
  code.list_size = 10;
  code.crc_dec = 7;
  code.crc_length = 3;
  code.generator_poly = {56721, 61713};

  std::vector<int> msg = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0};
  std::vector<int> crc_msg = calculateCRC(msg, code.crc_dec, code.crc_length);
  

  std::vector<int> checkMsg = {0, 0, 0, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 0, 0, 0};
  std::cout << "checking " << checkCRC(checkMsg, code.crc_dec, code.crc_length) << std::endl;

}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}