#include <gtest/gtest.h>
#include <vector>
#include "../src/feedForwardTrellis.h"
#include "../src/viterbiCodec.h"
#include "../src/minHeap.h"

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
  // std::vector<std::vector<Cell>> test_trellis = codec.getTrellis();

  // // stage 0
  // EXPECT_EQ(test_trellis[0][0].pathMetric,0);
  // // stage 1
  // EXPECT_EQ(test_trellis[0][1].pathMetric,2);
  // EXPECT_EQ(test_trellis[4][1].pathMetric,0);
  // EXPECT_EQ(test_trellis[6][1].pathMetric,INT_MAX);
  // // stage 2
  // EXPECT_EQ(test_trellis[0][2].pathMetric,3);
  // EXPECT_EQ(test_trellis[2][2].pathMetric,2);
  // EXPECT_EQ(test_trellis[4][2].pathMetric,3);
  // EXPECT_EQ(test_trellis[6][2].pathMetric,0);
  // // stage 3
  // EXPECT_EQ(test_trellis[0][3].pathMetric,4);
  // EXPECT_EQ(test_trellis[0][3].fatherState,0);
  // EXPECT_EQ(test_trellis[1][3].pathMetric,3);
  // EXPECT_EQ(test_trellis[1][3].fatherState,2);
  // EXPECT_EQ(test_trellis[3][3].pathMetric,0);
  // EXPECT_EQ(test_trellis[4][3].pathMetric,4);
  // // stage 4
  // EXPECT_EQ(test_trellis[0][4].pathMetric,3);
  // EXPECT_EQ(test_trellis[0][4].fatherState,1);
  // EXPECT_EQ(test_trellis[0][4].subPathMetric,6);
  // EXPECT_EQ(test_trellis[0][4].subFatherState,0);
  // EXPECT_EQ(test_trellis[1][4].pathMetric,2);
  // EXPECT_EQ(test_trellis[1][4].fatherState,3);
  // EXPECT_EQ(test_trellis[1][4].subPathMetric,5);
  // EXPECT_EQ(test_trellis[1][4].subFatherState,2);
  // EXPECT_EQ(test_trellis[4][4].pathMetric,4);
  // EXPECT_EQ(test_trellis[4][4].fatherState,0);
  // EXPECT_EQ(test_trellis[4][4].subPathMetric,5);
  // EXPECT_EQ(test_trellis[4][4].subFatherState,1);
  // // best path
  // EXPECT_EQ(test_trellis[5][4].pathMetric,0);
  // EXPECT_EQ(test_trellis[5][4].fatherState,3);
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
}


int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}