#include <gtest/gtest.h>
#include "../include/minHeap.h"
#include "../include/viterbiCodec.h"
#include "../include/feedForwardTrellis.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <random>

namespace {

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


template<typename T>
void outputMat(const std::vector<T>& vec, std::ofstream& outputFile) {
  outputFile << "[";
  for (const T& ele : vec) {
    outputFile << ele << ", ";
  }
  outputFile << "] with size: " << vec.size();
} 

}  // private namespace

TEST(ListDecoderTest, FloatPrecisionRanking) {

  std::ofstream outputFile("output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
  }
  std::vector<int> transmitted_msg = {1, 1, 1, 0, 1, 0, 1, 1};

  // modulated_signal
  std::vector<int> modulated_signal = {-1, -1, -1, 1, -1, 1, -1, -1};
  outputFile << "printing modulated signal: ";
  outputMat(modulated_signal, outputFile);
  outputFile << std::endl;

  double SNR = 10.0;
  // noisy signal
  std::vector<double> noisy_signal = addNoise(modulated_signal, SNR);
  
  outputFile << "printing noisy signal: ";
  outputMat(noisy_signal, outputFile);
  outputFile << std::endl;

  // DECODING
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 3;
  code.list_size = 10;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {13, 17};
  ViterbiCodec codec(code);
  std::vector<MessageInformation> output = codec.unconstraintListViterbiDecoding(noisy_signal);

  
  outputFile << "printing single 2^3 states list decoding results: " << std::endl;
  for (int i = 0; i < output.size(); ++i) {
    outputFile << " " << i << "th path  (" << "list rank = " << output[i].list_rank << "): ";
    outputMat(output[i].path, outputFile); 
    outputFile << " with pathMetric = " << output[i].path_metric;
    outputFile << std::endl;
  }

  outputFile.close();
}





int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}