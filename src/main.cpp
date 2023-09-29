#include <fstream>
#include <iostream>
#include <random>

#include "viterbiCodec.h"
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

} // namespace BPSK

int main() {
  // output path
  std::string outputFilePath = "../output/";
  std::ofstream outputFile(outputFilePath + "output.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return 1;
  }

  // SNR
  double SNR_dB_start = 0.0;
  double SNR_dB_end = 20.0;
  double SNR_dB_step = 2.0;
  std::vector<double> SNR_dB = {};
  std::vector<double> SNR = {};
  for (double i = SNR_dB_start; i <= SNR_dB_end; i += SNR_dB_step) {
    SNR_dB.push_back(i);
    SNR.push_back(pow(10.0, i / 10.0));
  }

  int mc_N = 1e4;
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 3;
  code.generator_poly = {13, 17};

  ViterbiCodec codec(code);
  
  int seed = 47;
  std::mt19937 generator(seed);

  int num_bits = 10; 
  std::vector<int> msg;
  for (int i = 0; i < num_bits; ++i) {
    int random_bit = generator() % 2;
    msg.push_back(random_bit);
  }
  outputFile << "Printing original message: " << std::endl;
  CodecUtils::output(msg, outputFile);

  // coding
  std::vector<int> encoded_msg = codec.encodeZTCC(msg);
  outputFile << "Printing coded message: " << std::endl;
  CodecUtils::output(encoded_msg, outputFile);

  std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
  outputFile << "Printing modulated signal: " << std::endl;
  CodecUtils::output(modulated_signal, outputFile);

  std::vector<double> received_signal = AWGN::addNoise(modulated_signal, 5.0);
  outputFile << "Printing received signal: " << std::endl;
  CodecUtils::output(received_signal, outputFile);

  outputFile << "measuring euclidean distance: " << CodecUtils::euclideanDistance(received_signal, modulated_signal)<< std::endl;


  outputFile.close();
  return 0;
}
