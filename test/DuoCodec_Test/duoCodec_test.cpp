#include <cassert>
#include <gtest/gtest.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <random>
#include "../include/feedForwardTrellis.h"
#include "../include/viterbiCodec.h"
#include "../include/minHeap.h"

namespace {

template<typename T>
void outputMat(const std::vector<T>& vec, std::ofstream& outputFile) {
  outputFile << "[";
  for (const T& ele : vec) {
    outputFile << ele << ", ";
  }
  outputFile << "] with size: " << vec.size();
} 


}  // private namespace
namespace AWGN {

std::vector<double> addNoise(std::vector<int> modulated_signal, double SNR) {
  std::random_device rd;
  std::mt19937 noise_gen( 47 );  // 82
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


}  // namespace AWGN

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


TEST(DuoCodecTest, EncodeTest) {
  
  std::ofstream outputFile("output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
  }
  CodeInformation code_1;
  // x^14+x^13+x^9+x^8+x^7+x^6+x^3+x+1 = 
  //                        CRC: (x^3+x^2+1) 
  //                        generator: (x^11+x^8+x^7+x^3+x^2+x+1)
  code_1.k = 1;
  code_1.n = 1;
  code_1.v = 11;
  code_1.list_size = 10;
  code_1.crc_dec = 13;
  code_1.crc_length = 4;
  code_1.generator_poly = {4617};

  CodeInformation code_2;
  // x^14+x^12+x^11+x^10+x^8+x^7+x^6+x^4+1 =
  //                         CRC: (x^2+x+1) 
  //                         generator: x^12+x^11+x^10+x^9+x^8+x^5+x^3+x+1 
  code_2.k = 1;
  code_2.n = 1;
  code_2.v = 12;
  code_2.list_size = 10;
  code_2.crc_dec = 7;
  code_2.crc_length = 3;
  code_2.generator_poly = {17453};

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  CodeInformation code_3;
  code_3.k = 1;
  code_3.n = 2;
  code_3.v = 14;
  code_3.list_size = 100;
  code_3.crc_dec = -1;
  code_3.crc_length = -1;
  code_3.generator_poly = {56721, 61713};
  ViterbiCodec codec_3(code_3);


  double snr = 2.0;

  // ZTCC Encode Test
  outputFile << "Running simulation under snr = " << snr << std::endl;

  // generate random message
  std::vector<int> msg = {1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1};
  outputFile << "msg: ";
  outputMat(msg, outputFile);
  outputFile << std::endl;


  // encode crc convolved message using a rate (1,2) convolutional code
  std::vector<int> encoded_msg = codec_3.encodeZTCC(msg);
  outputFile << "encoded_msg: ";
  outputMat(encoded_msg, outputFile);
  outputFile << std::endl;

  assert(encoded_msg.size() == (msg.size() + code_3.v) * 2);
  
  // code 1
  std::vector<int> crc_msg_1 = codec_1.calculateCRC(msg);
  outputFile << "crc_msg_1: ";
  outputMat(crc_msg_1, outputFile);
  outputFile << std::endl;

  std::vector<int> crc_msg_11 = codec_1.convolveCRC(msg);
  outputFile << "crc_msg_11: ";
  outputMat(crc_msg_11, outputFile);
  outputFile << std::endl;

  std::vector<int> encoded_crc_msg_1 = codec_1.encodeZTCC(crc_msg_1);
  outputFile << "encoded_crc_msg_1: ";
  outputMat(encoded_crc_msg_1, outputFile);
  outputFile << std::endl;

  std::vector<int> encoded_crc_msg_11 = codec_1.encodeZTCC(crc_msg_11);
  outputFile << "encoded_crc_msg_11: ";
  outputMat(encoded_crc_msg_11, outputFile);
  outputFile << std::endl;

  std::vector<int> modulated_signal = BPSK::modulate(encoded_crc_msg_11);
  outputFile << "modulated signal 1: ";
  outputMat(modulated_signal, outputFile);
  outputFile << std::endl;

  std::vector<double> received_signal = AWGN::addNoise(modulated_signal, snr);
  outputFile << "Received signal 1: ";
  outputMat(received_signal, outputFile);
  outputFile << std::endl;
  
  outputFile << std::endl;

  // code 2
  std::vector<int> crc_msg_2 = codec_2.calculateCRC(msg);
  outputFile << "crc_msg_2: ";
  outputMat(crc_msg_2, outputFile);
  outputFile << std::endl;

  std::vector<int> crc_msg_22 = codec_2.convolveCRC(msg);
  outputFile << "crc_msg_22: ";
  outputMat(crc_msg_22, outputFile);
  outputFile << std::endl;

  std::vector<int> encoded_crc_msg_2 = codec_2.encodeZTCC(crc_msg_2);
  outputFile << "encoded_crc_msg_2: ";
  outputMat(encoded_crc_msg_2, outputFile);
  outputFile << std::endl;

  std::vector<int> encoded_crc_msg_22 = codec_2.encodeZTCC(crc_msg_22);
  outputFile << "encoded_crc_msg_22: ";
  outputMat(encoded_crc_msg_22, outputFile);
  outputFile << std::endl;

  std::vector<int> modulated_signal_2 = BPSK::modulate(encoded_crc_msg_22);
  outputFile << "modulated signal 2: ";
  outputMat(modulated_signal_2, outputFile);
  outputFile << std::endl;

  std::vector<double> received_signal_2 = AWGN::addNoise(modulated_signal_2, snr);
  outputFile << "Received signal 2: ";
  outputMat(received_signal_2, outputFile);
  outputFile << std::endl;
  outputFile << std::endl;

  // code 3
  std::vector<int> encoded_msg_14 = codec_3.encodeZTCC(msg);
  outputFile << "encoded_msg_14: ";
  outputMat(encoded_msg_14, outputFile);
  outputFile << std::endl;

  std::vector<int> modulated_signal_14 = BPSK::modulate(encoded_msg_14);
  outputFile << "modulated signal 14: ";
  outputMat(modulated_signal_14, outputFile);
  outputFile << std::endl;

  std::vector<double> received_signal_14 = AWGN::addNoise(modulated_signal_14, snr);
  outputFile << "Received signal 14: ";
  outputMat(received_signal_14, outputFile);
  outputFile << std::endl;

  std::vector<double> received_codec_2;
  std::vector<double> received_codec_1;

  // unleaver to unleave the bits from received_signal_14
  for (size_t i = 0; i < received_signal_14.size(); ++i) {
        if (i % 2 == 0) {
            received_codec_2.push_back(received_signal_14[i]);
        } else {
            received_codec_1.push_back(received_signal_14[i]);
        }
  }
  
  outputFile << "For codec 1: ";
  outputMat(received_codec_1, outputFile);
  outputFile << std::endl;

  outputFile << "For codec 2: ";
  outputMat(received_codec_2, outputFile);
  outputFile << std::endl;

  // Decoding starts
  outputFile << "Decoding begins -----------------" << std::endl;

  std::vector<MessageInformation> output_1 = codec_1.ZTCCListViterbiDecoding(received_codec_1);
  std::vector<MessageInformation> output_2 = codec_2.ZTCCListViterbiDecoding(received_codec_2);

  std::vector<MessageInformation> output_3 = codec_3.unconstraintZTCCDecoding(received_signal_14);

  outputFile << std::endl;

  outputFile << "Printing Duo list decoder results: " << std::endl;
  for (int i = 0; i < output_1.size(); ++i) {
    outputFile << " " << i << "th message  (" << "list rank = " << output_1[i].list_rank << "): ";
    outputMat(output_1[i].message, outputFile); 
    outputFile << " with pathMetric = " << output_1[i].path_metric;
    outputFile <<  "     " << "(list rank = " << output_2[i].list_rank << "): ";
    outputMat(output_2[i].message, outputFile);
    outputFile << " with pathMetric = " << output_2[i].path_metric;
    outputFile << std::endl;
  }

  outputFile << std::endl;

  outputFile << "printing single 2^14 states list decoding results: " << std::endl;
  for (int i = 0; i < output_3.size(); ++i) {
    outputFile << " " << i << "th message  (" << "list rank = " << output_3[i].list_rank << "): ";
    outputMat(output_3[i].message, outputFile); 
    outputFile << " with pathMetric = " << output_3[i].path_metric;
    outputFile << std::endl;
  }




  outputFile.close();
}

// TEST(DuoCodecTest, TBCCDecodingCodec1) {
//   std::ofstream outputFile("output.txt");
//   if (!outputFile.is_open()) {
//     std::cerr << "Failed to open the file for writing." << std::endl;
//   }
//   CodeInformation code_1;
//   // x^14+x^13+x^9+x^8+x^7+x^6+x^3+x+1 = 
//   //                        CRC: (x^3+x^2+1) 
//   //                        generator: (x^11+x^8+x^7+x^3+x^2+x+1)
//   code_1.k = 1;
//   code_1.n = 1;
//   code_1.v = 11;
//   code_1.list_size = 10;
//   code_1.crc_dec = 13;
//   code_1.crc_length = 4;
//   code_1.generator_poly = {4617};

//   CodeInformation code_2;
//   // x^14+x^12+x^11+x^10+x^8+x^7+x^6+x^4+1 =
//   //                         CRC: (x^2+x+1) 
//   //                         generator: x^12+x^11+x^10+x^9+x^8+x^5+x^3+x+1 
//   code_2.k = 1;
//   code_2.n = 1;
//   code_2.v = 12;
//   code_2.list_size = 10;
//   code_2.crc_dec = 7;
//   code_2.crc_length = 3;
//   code_2.generator_poly = {17453};

//   ViterbiCodec codec_1(code_1);
//   ViterbiCodec codec_2(code_2);



// }




int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}