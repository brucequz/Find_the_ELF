#include <fstream>
#include <iostream>
#include <random>

#include "../include/viterbiCodec.h"
#include "mat.h"
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

namespace MatlabUtils {

template <typename T>
void writeToMat(const std::vector<std::vector<T>>& data, const char* filePath, const char* varName) {
  // Open MAT file for writing
  MATFile* pmat = matOpen(filePath, "w");
  if (pmat == NULL) {
      std::cerr << "Error creating MAT file." << std::endl;
      return;
  }

  // Get the number of rows and columns in the data
  mwSize rows = static_cast<mwSize>(data.size());
  mwSize cols = (rows > 0) ? static_cast<mwSize>(data[0].size()) : 0;

  // Create a mxArray and populate it with data
  mxArray* pmxArray = mxCreateDoubleMatrix(rows, cols, mxREAL);
  if (pmxArray == NULL) {
      std::cerr << "Error creating mxArray." << std::endl;
      matClose(pmat);
      return;
  }

  // Populate the mxArray with data
  double* pr = mxGetPr(pmxArray);
  for (mwIndex i = 0; i < rows; ++i) {
      for (mwIndex j = 0; j < cols; ++j) {
          pr[i + j * rows] = data[i][j];
      }
  }

  // Write the mxArray to the MAT file with the specified variable name
  int status = matPutVariable(pmat, varName, pmxArray);
  if (status != 0) {
      std::cerr << "Error writing mxArray to MAT file." << std::endl;
      mxDestroyArray(pmxArray);
      matClose(pmat);
      return;
  }

  // Close MAT file and clean up resources
  mxDestroyArray(pmxArray);
  matClose(pmat);
}

} // namespace MatlabUtils

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

  int mc_N = 2;
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
  
  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 15; 
  for (int i = 0; i < mc_N; ++i) {
    for (double snr_dB : {0.0}) {
      
      outputFile << "Now working on snr: " << snr_dB << "-------------------" << std::endl;

      std::vector<int> msg;
      for (int i = 0; i < num_bits; ++i) {
        int random_bit = msg_gen() % 2;
        msg.push_back(random_bit);
      }
      outputFile << "Printing original message: " << std::endl;
      CodecUtils::output(msg, outputFile);
      outputFile << std::endl;

      // add crc
      std::vector<int> crc_msg = codec.calculateCRC(msg);
      outputFile << crc_msg.size() << " Printing crc message: " << std::endl;
      CodecUtils::output(crc_msg, outputFile);
      outputFile << std::endl;

      // coding
      std::vector<int> encoded_msg = codec.encodeZTCC(crc_msg);
      outputFile << encoded_msg.size() << " Printing coded message: " << std::endl;
      CodecUtils::output(encoded_msg, outputFile);
      outputFile << std::endl;

      std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
      outputFile << "Printing modulated signal: " << std::endl;
      CodecUtils::output(modulated_signal, outputFile);
      outputFile << std::endl;

      std::vector<double> received_signal = AWGN::addNoise(modulated_signal, snr_dB);
      outputFile << "Printing received signal: " << std::endl;
      CodecUtils::output(received_signal, outputFile);
      outputFile << std::endl;

      outputFile << "measuring euclidean distance: " << CodecUtils::euclideanDistance(received_signal, modulated_signal)<< std::endl;

      std::vector<int> demodulated_signal = BPSK::demodulate(received_signal);
      outputFile << "Printing demodulated signal: " << std::endl;
      CodecUtils::output(demodulated_signal, outputFile);
      outputFile << std::endl;

      std::vector<int> decoded_msg = codec.softViterbiDecode(received_signal).message;
      outputFile << "Printing soft decoded message: " << std::endl;
      CodecUtils::output(decoded_msg, outputFile);
      outputFile << std::endl;

      std::vector<int> hard_decoded_msg = codec.viterbiDecode(demodulated_signal).message;
      outputFile << "Printing hard decoded message: " << std::endl;
      CodecUtils::output(hard_decoded_msg, outputFile);
      outputFile << std::endl;

      std::vector<MessageInformation> output = codec.listViterbiDecoding(received_signal);
      outputFile << "Printing list decoder results: " << std::endl;
      for (int i = 0; i < output.size(); ++i) {
        outputFile << output[i].message.size() << " " << i << "th message:  ";
        CodecUtils::output(output[i].message, outputFile);
        outputFile << std::endl;
      }
      
      // truncate the flushing bits
      decoded_msg.resize(msg.size());
      hard_decoded_msg.resize(msg.size());
      outputFile << "Measuring hamming distance: " << CodecUtils::hammingDistance(msg, hard_decoded_msg) << std::endl;

    }
  }

  outputFile.close();
  return 0;
}
