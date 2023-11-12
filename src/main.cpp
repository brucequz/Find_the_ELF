#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>

#include "../include/dualListDecoder.h"
#include "../include/viterbiCodec.h"
// #include "mat.h"

static std::random_device rd{};
static std::mt19937 noise_gen(rd());  // 82

namespace AWGN {

std::vector<double> addNoise(std::vector<int> modulated_signal, double SNR) {
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

// namespace MatlabUtils {

// template <typename T>
// void writeToMat(const std::vector<std::vector<T>>& data, const char* filePath,
//                 const char* varName) {
//   // Open MAT file for writing
//   MATFile* pmat = matOpen(filePath, "w");
//   if (pmat == NULL) {
//     std::cerr << "Error creating MAT file." << std::endl;
//     return;
//   }

//   // Get the number of rows and columns in the data
//   mwSize rows = static_cast<mwSize>(data.size());
//   mwSize cols = (rows > 0) ? static_cast<mwSize>(data[0].size()) : 0;

//   // Create a mxArray and populate it with data
//   mxArray* pmxArray = mxCreateDoubleMatrix(rows, cols, mxREAL);
//   if (pmxArray == NULL) {
//     std::cerr << "Error creating mxArray." << std::endl;
//     matClose(pmat);
//     return;
//   }

//   // Populate the mxArray with data
//   double* pr = mxGetPr(pmxArray);
//   for (mwIndex i = 0; i < rows; ++i) {
//     for (mwIndex j = 0; j < cols; ++j) {
//       pr[i + j * rows] = data[i][j];
//     }
//   }

//   // Write the mxArray to the MAT file with the specified variable name
//   int status = matPutVariable(pmat, varName, pmxArray);
//   if (status != 0) {
//     std::cerr << "Error writing mxArray to MAT file." << std::endl;
//     mxDestroyArray(pmxArray);
//     matClose(pmat);
//     return;
//   }

//   // Close MAT file and clean up resources
//   mxDestroyArray(pmxArray);
//   matClose(pmat);
// }

// }  // namespace MatlabUtils

int main() {
  // output path
  std::string outputFilePath = "../output/even_smaller_example/";
  std::ofstream outputFile(outputFilePath + "std_snr_10-50.txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
    return 1;
  }

  std::vector<double> EbN0 = {1, 2, 2.5, 3, 3.5, 4, 4.5, 5};
  std::vector<double> SNR_dB;  // SNR is required for noise computations
  double offset = 10 * log10((double)2 * 64 /
                             (double)134);  // real rate of this code is 32/512
  for (int i = 0; i < EbN0.size(); i++) {
    SNR_dB.push_back(EbN0[i]);
  }

  // SNR
  // double SNR_dB_start = 1.0;
  // double SNR_dB_end = 7.0;
  // double SNR_dB_step = 0.5;
  // std::vector<double> SNR_dB = {};
  // std::vector<double> SNR = {};
  // for (double i = SNR_dB_start; i <= SNR_dB_end; i += SNR_dB_step) {
  //   SNR_dB.push_back(i);
  //   SNR.push_back(pow(10.0, i / 10.0));
  // }

  // int mc_N = 10000;
  int max_errors = 1000;
  int list_size = 1;

  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 3;
  code.list_size = 1;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {13, 17};  // octal
  ViterbiCodec codec(code);

  CodeInformation code_1;
  // 61713
  // x^14+x^13+x^9+x^8+x^7+x^6+x^3+x+1 =
  //                        CRC: (x^3+x^2+1)
  //                        generator: (x^11+x^8+x^7+x^3+x^2+x+1)

  // 3217
  // x^10 + x^9 + x^7 + x^3 + x^2 + x + 1
  //                        CRC: (x^4+x^2+1)
  //                        generator: x^6+x^5+x^4+x+1
  code_1.k = 1;
  code_1.n = 1;
  code_1.v = 6;
  code_1.list_size = list_size;
  code_1.crc_dec = 21;
  code_1.crc_length = 5;
  code_1.generator_poly = {163};  // octal

  CodeInformation code_2;
  // 56721
  // x^14+x^12+x^11+x^10+x^8+x^7+x^6+x^4+1 =
  //                         CRC: (x^2+x+1)
  //                         generator: x^12+x^11+x^10+x^9+x^8+x^5+x^3+x+1

  // 2473
  // x^10 + x^8 + x^5 + x^4 + x^3 + x + 1
  //                         CRC: (x^3+x+1)
  //                         generator: x^7+x^4+1
  code_2.k = 1;
  code_2.n = 1;
  code_2.v = 7;
  code_2.list_size = list_size;
  code_2.crc_dec = 11;
  code_2.crc_length = 4;
  code_2.generator_poly = {221};

  //ViterbiCodec codec_1(code_1);
  //ViterbiCodec codec_2(code_2);

  // std::vector<CodeInformation> dld_codes = {code_1, code_2};
  // DualListDecoder DLD(dld_codes, 1000000);

  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 64;

  std::vector<double> correct_decoding_snr;
  std::vector<double> ML_decoding_error_snr;

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_error_vec;

  // matlab debug purposes
  // std::vector<std::vector<int>> matlab_message;
  // std::vector<std::vector<double>> matlab_received_signal;

  for (double snr_dB : SNR_dB) {
    std::cout << "Now working on snr: " << snr_dB << "-------------------"
              << std::endl;

    std::vector<int> expected_list_ranks = {0, 0};

    // standard viterbi
    int ML_decoding_error = 0;
    int Correct_decoding = 0;

    // DLD decoding
    int DLD_correct = 0;
    int DLD_error = 0;

    int number_of_trials = 0;
    int number_of_errors = 0;

    while (number_of_errors < max_errors) {
      number_of_trials++;

      if (number_of_trials % 5000 == 0) {
        std::cout << "Trial number: " << number_of_trials << std::endl;
        std::cout << "Current number of errors: " << number_of_errors
                  << std::endl;
      }

      std::vector<int> msg;
      for (int i = 0; i < num_bits; ++i) {
        int random_bit = rand() % 2;
        msg.push_back(random_bit);
      }

      // coding
      std::vector<int> encoded_msg = codec.encodeZTCC(msg);
      outputFile << " Printing coded message: " << std::endl;
      CodecUtils::outputMat(encoded_msg, outputFile);
      outputFile << std::endl;

      assert(encoded_msg.size() == (msg.size() + code.v) * 2);

      std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
      // outputFile << "Printing modulated signal: " << std::endl;
      // CodecUtils::outputMat(modulated_signal, outputFile);
      // outputFile << std::endl;

      std::vector<double> received_signal =
          AWGN::addNoise(modulated_signal, snr_dB);
      outputFile << "Printing received signal: " << std::endl;
      CodecUtils::outputMat(received_signal, outputFile);
      outputFile << std::endl;

      // Creating matlab file for a single snr
      // matlab_message.push_back(msg);
      // matlab_received_signal.push_back(received_signal);

      // SOFT VITERBI DECODING
      MessageInformation output_ML = codec.softViterbiDecodeZTCC(received_signal);

      assert(output_ML.message.size() == msg.size());
      if (CodecUtils::areVectorsEqual(output_ML.message, msg)) {
        Correct_decoding++;
      } else {
        ML_decoding_error++;
        number_of_errors++;
      }

      // // DLD DECODING
      // DLDInfo output_DLD = DLD.adaptiveDecode(received_signal);
      // if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
      //   DLD_correct++;
      // } else {
      //   DLD_error++;
      // }

      // // update list ranks
      // for (int i = 0; i < expected_list_ranks.size(); ++i) {
      //   expected_list_ranks[i] += output_DLD.list_ranks[i];
      // }

      // update in the outputfile
      if (number_of_trials % 2000 == 0) {
        outputFile << "Trial: " << number_of_trials << std::endl;
      }
    }

    // for (int j = 0; j < expected_list_ranks.size(); ++j) {
    //   expected_list_ranks[j] /= mc_N;
    // }
    std::cout << "For snr = " << snr_dB << ", " << number_of_trials
              << " were ran to accumulate " << number_of_errors << " errors."
              << std::endl;
    std::cout << "std viterbi correct decoding: " << Correct_decoding
              << " , Percentage: "
              << (double)Correct_decoding / number_of_trials << std::endl;
    std::cout << "std viterbi wrong decoding: " << ML_decoding_error
              << " , Percentage: "
              << (double)ML_decoding_error / number_of_trials << std::endl;
    // std::cout << "DLD correct decoding: " << DLD_correct << " , Percentage: "
    // << (double)DLD_correct/mc_N << std::endl; std::cout << "DLD wrong
    // decoding: " << DLD_error << " , Percentage: " << (double)DLD_error/mc_N
    // << std::endl;

    // std::cout << "expected list ranks: [ " << expected_list_ranks[0] << ", "
    // << expected_list_ranks[1] << " ]" << std::endl;

    correct_decoding_snr.push_back((double)Correct_decoding / number_of_trials);
    ML_decoding_error_snr.push_back((double)ML_decoding_error /
                                    number_of_trials);
    // DLD_correct_vec.push_back((double)DLD_correct/mc_N);
    // DLD_error_vec.push_back((double)DLD_error/mc_N);
  }
  // MatlabUtils::writeToMat(matlab_message, "../test/Matlab_decoding_acc_verification/messages.mat", "messages");
  // MatlabUtils::writeToMat(matlab_received_signal, "../test/Matlab_decoding_acc_verification/received_signal.mat", "received_signal");

  std::cout << "Std viterbi CORRECT decoding percentage: ";
  CodecUtils::print(correct_decoding_snr);
  std::cout << std::endl;
  std::cout << "Std viterbi WRONG decoding percentage: ";
  CodecUtils::print(ML_decoding_error_snr);
  std::cout << std::endl;
  // std::cout << "DLD Decoder CORRECT decoding percentage: ";
  // CodecUtils::print(DLD_correct_vec); std::cout << std::endl; std::cout <<
  // "DLD Decoder WRONG decoding percentage: ";
  // CodecUtils::print(DLD_error_vec); std::cout << std::endl;
  outputFile.close();
  return 0;
}
