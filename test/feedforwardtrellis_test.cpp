#include <gtest/gtest.h>
#include <vector>
#include <random>
#include "../src/feedForwardTrellis.h"
#include "../src/viterbiCodec.h"
#include "../src/minHeap.h"

namespace AWGN {

std::vector<double> addNoise(std::vector<int> modulated_signal, double SNR) {
  std::mt19937 generator;
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
    noisyMsg.push_back(modulated_signal[i] + distribution(generator));
  }
  return noisyMsg;
}

}  // namespace AWGN

namespace BPSK {
  
std::vector<int> modulate(std::vector<int> encoded_msg) {
  std::vector<int> modulated_signal(encoded_msg.size());
  for (int i = 0; i < encoded_msg.size(); ++i) {
    modulated_signal[i] = -2 * encoded_msg[i] + 1;
  }
  return modulated_signal;
}

}  // namespace BPSK

namespace Utils {

template <typename T>
void print(const std::vector<T>& vec) {
  for (const T& element : vec) {
      std::cout << element << " ";
  }
  std::cout << std::endl;
}

template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
  for (const std::vector<T>& row : matrix) {
      for (const T& element : row) {
          std::cout << element << " ";
      }
      std::cout << std::endl;
  }
}

}

// Define test cases
TEST(TrellisTest, DecodeNoError) {
  std::vector<int> poly = {13, 17};
  FeedForwardTrellis trellis(1, 2, 3, poly);
  std::vector<int> input_1 = {1, 0, 1, 0, 1, 1, 1};
  std::vector<int> encoded = trellis.encode(input_1);
  std::vector<int> expect_encoded = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 1};
  std::vector<int> input_2 = {0, 0, 1, 0, 1, 1, 0};
  std::vector<int> encoded_2 = trellis.encode(input_2);
  std::vector<int> expect_encoded_2 = {0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 0};
  std::vector<int> nonexpect_encoded_3 = {0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 1, 1, 1};
  
  EXPECT_EQ(encoded,expect_encoded);
  EXPECT_EQ(encoded_2,expect_encoded_2);
  EXPECT_NE(encoded_2,nonexpect_encoded_3);
}

TEST(CodecTest, BuildDecodeTrellis) {
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);
  std::vector<int> received_message = {1, 1, 1, 0, 1, 0, 1, 1};
  codec.viterbiDecode(received_message);
  
}

TEST(CodecTest, DecodingPath) {
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);
  std::vector<int> received_message = {1, 1, 1, 0, 1, 0, 1, 1};
  MessageInformation output = codec.viterbiDecode(received_message);

  EXPECT_EQ(output.path.back(), 5);
  EXPECT_EQ(output.path[3], 3);
  EXPECT_EQ(output.path[2], 6);
  EXPECT_EQ(output.path[1], 4);
  EXPECT_EQ(output.path[0], 0);
  EXPECT_EQ(output.message[0], 1);
  EXPECT_EQ(output.message[1], 1);
  EXPECT_EQ(output.message[2], 0);
  EXPECT_EQ(output.message[3], 1);
}

TEST(CodecTest, RandomDecoding) {
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);
  std::vector<int> received_message = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0};
  MessageInformation output = codec.viterbiDecode(received_message);

  // EXPECT_EQ(output.path.back(), 5);
  // EXPECT_EQ(output.path[3], 3);
  // EXPECT_EQ(output.path[2], 6);
  // EXPECT_EQ(output.path[1], 4);
  // EXPECT_EQ(output.path[0], 0);
  EXPECT_EQ(output.message[0], 1);
  EXPECT_EQ(output.message[1], 0);
  EXPECT_EQ(output.message[2], 1);
  EXPECT_EQ(output.message[3], 0);
  EXPECT_EQ(output.message[4], 1);

  std::vector<int> received_message2 = {1, 1, 1, 0, 0, 1, 1, 0, 1, 0};
  MessageInformation output2 = codec.viterbiDecode(received_message2);
  EXPECT_EQ(output2.message[0], 1);
  EXPECT_EQ(output2.message[1], 1);
  EXPECT_EQ(output2.message[2], 1);
  EXPECT_EQ(output2.message[3], 1);
  EXPECT_EQ(output2.message[4], 1);
}

TEST(CodecTest, LargeMemoryElements) {
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);
  // std::vector<int> poly = {56721, 61713};
  // ViterbiCodec codec(1, 2, 14, poly);
  std::vector<int> received_message = {1, 1, 0, 1, 0, 0, 1, 0, 0, 0};
  MessageInformation output = codec.viterbiDecode(received_message);

  // EXPECT_EQ(output.path.back(), 5);
  // EXPECT_EQ(output.path[3], 3);
  // EXPECT_EQ(output.path[2], 6);
  // EXPECT_EQ(output.path[1], 4);
  // EXPECT_EQ(output.path[0], 0);
  EXPECT_EQ(output.message[0], 1);
  EXPECT_EQ(output.message[1], 0);
  EXPECT_EQ(output.message[2], 1);
  EXPECT_EQ(output.message[3], 0);
  EXPECT_EQ(output.message[4], 1);

  std::vector<int> received_message2 = {1, 1, 1, 0, 0, 1, 1, 0, 1, 0};
  MessageInformation output2 = codec.viterbiDecode(received_message2);
  EXPECT_EQ(output2.message[0], 1);
  EXPECT_EQ(output2.message[1], 1);
  EXPECT_EQ(output2.message[2], 1);
  EXPECT_EQ(output2.message[3], 1);
  EXPECT_EQ(output2.message[4], 1);
}

TEST(CodecTest, ZTCCEncodeTest) {
  std::vector<int> poly = {13, 17};
  ViterbiCodec codec(1, 2, 3, poly);
  std::vector<int> message = {1, 0, 1, 1};
  std::vector<int> encoded = codec.encodeZTCC(message);

  EXPECT_EQ(encoded.size(), 14);
  EXPECT_EQ(encoded[0], 1);
  EXPECT_EQ(encoded[1], 1);
  EXPECT_EQ(encoded[2], 0);
  EXPECT_EQ(encoded[3], 1);
  EXPECT_EQ(encoded[4], 0);
  EXPECT_EQ(encoded[5], 0);
  EXPECT_EQ(encoded[6], 0);
  EXPECT_EQ(encoded[7], 1);
  EXPECT_EQ(encoded[8], 1);
  EXPECT_EQ(encoded[9], 0);
  EXPECT_EQ(encoded[10], 0);
  EXPECT_EQ(encoded[11], 0);
  EXPECT_EQ(encoded[12], 1);
  EXPECT_EQ(encoded[13], 1);
  
  std::vector<int> modulated_signal = BPSK::modulate(encoded);
  std::cout << "printing modualted signal: " << std::endl;
  Utils::print(modulated_signal);

  std::vector<double> noisy_signal = AWGN::addNoise(modulated_signal, 0.0);
  std::cout << "printing noisy signal: " << std::endl;
  Utils::print(noisy_signal);
  
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}