#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <map>
#include <queue>
#include <random>

#include "../include/dualListDecoder.h"
#include "../include/viterbiCodec.h"
// #include "mat.h"

static std::random_device rd{};
static std::mt19937 noise_gen(82);  // 82

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

/*
Functions that end with Experiment is a stand-alone experiment designed to run by itself
Functions taht end with Decoding is a decoding function that takes in a single received message and 
*/
void dualListExperiment(std::vector<double> SNR_dB, int maximum_list_size, int max_errors);
void softViterbiExperiment(std::vector<double> SNR_dB, int max_errors);



int main(int argc, char* argv[]) {

  // variables setup
  // 1. snr
  std::vector<double> EbN0;
  std::istringstream iss(argv[1]);
  std::string token;

  while (std::getline(iss, token, ',')) {
      // Convert each substring to a double and add it to the vector
      EbN0.push_back(std::stod(token));
  }
  // std::vector<double> EbN0 = {3, 3.5};
  std::vector<double> SNR_dB;  // SNR is required for noise computations
  double offset = 10 * log10((double)2 * 64 /
                             (double)134);  // real rate of this code is 32/512
  for (int i = 0; i < EbN0.size(); i++) {
    SNR_dB.push_back(EbN0[i]);
  }

  // 2. number of maximum errors to accumulate
  int max_errors = 50;
  int list_size = 1;
  

  // 3. DLD decoding setup
  int DLD_maximum_list_size = 1e6;
  dualListExperiment(SNR_dB, DLD_maximum_list_size, max_errors);

  // soft Viterbi Experiment Setup
  // softViterbiExperiment(SNR_dB, max_errors);
  return 0;
}

void softViterbiExperiment(std::vector<double> SNR_dB, int max_errors) {
  // output path
  std::string outputFilePath = "../output/smaller_example/";
  std::ofstream outputFile(outputFilePath + "STD_snr_" + std::to_string(int(SNR_dB.front())) + "-" + 
                           std::to_string(int(SNR_dB.back())) + ".txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
  }

  /// Building Codec
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 10;
  code.list_size = 1;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {2473, 3217};  // octal
  ViterbiCodec codec(code);
  
  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 64;

  std::vector<double> correct_decoding_snr;
  std::vector<double> ML_decoding_error_snr;
  
  std::cout << "Running STD tests for v = " << code.v << " test, with generator poly: " << code.generator_poly[0] << ", " << 
  code.generator_poly[1] << " SNR = " << SNR_dB[0] << " - " << SNR_dB.back() << std::endl;

  for (double snr_dB : SNR_dB) {
    std::cout << "Now working on snr: " << snr_dB << "-------------------"
              << std::endl;

    std::vector<int> expected_list_ranks = {0, 0};
    
    int ML_decoding_error = 0;
    int Correct_decoding = 0;

    int number_of_trials = 0;
    int number_of_errors = 0;

    while (number_of_errors < max_errors) {
      number_of_trials++;

      if (number_of_trials % 500 == 0) {
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
      // outputFile << " Printing coded message: " << std::endl;
      // CodecUtils::outputMat(encoded_msg, outputFile);
      // outputFile << std::endl;

      assert(encoded_msg.size() == (msg.size() + code.v) * 2);

      std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
      // outputFile << "Printing modulated signal: " << std::endl;
      // CodecUtils::outputMat(modulated_signal, outputFile);
      // outputFile << std::endl;

      std::vector<double> received_signal =
          AWGN::addNoise(modulated_signal, snr_dB);
      // outputFile << "Printing received signal: " << std::endl;
      // CodecUtils::outputMat(received_signal, outputFile);
      // outputFile << std::endl;

      // Creating matlab file for a single snr
      // matlab_message.push_back(msg);
      // matlab_received_signal.push_back(received_signal);

      // // SOFT VITERBI DECODING
      MessageInformation output_ML = codec.softViterbiDecoding(received_signal);

      assert(output_ML.message.size() == msg.size());
      if (CodecUtils::areVectorsEqual(output_ML.message, msg)) {
        Correct_decoding++;
      } else {
        ML_decoding_error++;
        number_of_errors++;
      }

      // update in the outputfile
      if (number_of_trials % 2000 == 0) {
        outputFile << "Trial: " << number_of_trials << std::endl;
      }
    }

    std::cout << "For snr = " << snr_dB << ", " << number_of_trials
              << " were ran to accumulate " << number_of_errors << " errors."
              << std::endl;
    std::cout << "std viterbi correct decoding: " << Correct_decoding
              << " , Percentage: "
              << (double)Correct_decoding / number_of_trials << std::endl;
    std::cout << "std viterbi wrong decoding: " << ML_decoding_error
              << " , Percentage: "
              << (double)ML_decoding_error / number_of_trials << std::endl;
    correct_decoding_snr.push_back((double)Correct_decoding / number_of_trials);
    ML_decoding_error_snr.push_back((double)ML_decoding_error /
                                    number_of_trials);
  }


  ///////////  End of experiment ////////////
  // Matlab debug purposes
  // MatlabUtils::writeToMat(matlab_message, "../test/Matlab_decoding_acc_verification/messages.mat", "messages");
  // MatlabUtils::writeToMat(matlab_received_signal, "../test/Matlab_decoding_acc_verification/received_signal.mat", "received_signal");

  std::cout << "Std viterbi CORRECT decoding percentage: ";
  CodecUtils::print(correct_decoding_snr);
  std::cout << std::endl;
  std::cout << "Std viterbi WRONG decoding percentage: ";
  CodecUtils::print(ML_decoding_error_snr);
  std::cout << std::endl;
  outputFile.close();
}

void dualListExperiment(std::vector<double> SNR_dB, int maximum_list_size, int max_errors) {
  /*
  function description here
  TODO:
  
  */
  // output path
  std::string outputFilePath = "../output/smaller_example/";
  std::ostringstream oss;
  oss << std::scientific << maximum_list_size;

  std::string scientificString = oss.str();
  std::ofstream outputFile(outputFilePath + "DLD_" + scientificString + "_" + std::to_string(int(SNR_dB.front())) + "-" + 
                           std::to_string(int(SNR_dB.back())) + ".txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
  }
  

  /// Building Codec
  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 10;
  code.list_size = 1;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {2473, 3217};  // octal
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
  code_1.list_size = 1;
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
  code_2.list_size = 1;
  code_2.crc_dec = 11;
  code_2.crc_length = 4;
  code_2.generator_poly = {221};

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  std::vector<CodeInformation> dld_codes = {code_1, code_2};
  DualListDecoder DLD(dld_codes, maximum_list_size);


  // Simulation begin
  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 64;

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_list_exceed_vec;
  std::vector<double> DLD_error_vec;

  std::cout << "Running DLD tests for v = " << code.v << " test, with generator poly: " << code.generator_poly[0] << ", " << 
  code.generator_poly[1] << " SNR = " << SNR_dB[0] << " - " << SNR_dB.back() << std::endl;


  for (double snr_dB : SNR_dB) {
    std::cout << "Now working on snr: " << snr_dB << "-------------------"
              << std::endl;

    std::vector<int> expected_list_ranks = {0, 0};

    int DLD_correct = 0;
    int DLD_list_exceeded = 0;
    int DLD_error = 0;

    // vector to keep track of list sizes for both decoders
    std::vector<int> DLD_list_0_size;
    std::vector<int> DLD_list_1_size;

    int number_of_trials = 0;
    int number_of_errors = 0;

    while (number_of_errors < max_errors) {
      number_of_trials++;

      if (number_of_trials % 20 == 0) {
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
      // outputFile << " Printing coded message: " << std::endl;
      // CodecUtils::outputMat(encoded_msg, outputFile);
      // outputFile << std::endl;

      assert(encoded_msg.size() == (msg.size() + code.v) * 2);

      std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
      // outputFile << "Printing modulated signal: " << std::endl;
      // CodecUtils::outputMat(modulated_signal, outputFile);
      // outputFile << std::endl;

      std::vector<double> received_signal =
          AWGN::addNoise(modulated_signal, snr_dB);
      // outputFile << "Printing received signal: " << std::endl;
      // CodecUtils::outputMat(received_signal, outputFile);
      // outputFile << std::endl;

      // DLD DECODING
      DLDInfo output_DLD = DLD.adaptiveDecode(received_signal);
      // we save the list ranks
      DLD_list_0_size.push_back(output_DLD.list_ranks[0]);
      DLD_list_1_size.push_back(output_DLD.list_ranks[1]);

      if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
        // correct decoding
        DLD_correct++;
      } else if (output_DLD.message == std::vector<int>(64, -1)) {
        // list size exceeded
        // we save the message, received_signal and record the maximum list size
        DLD_list_exceeded++;
        number_of_errors++;
        outputFile << "LSE Error: " << std::endl;
        outputFile << "  Transmitted message: ";
        CodecUtils::outputMat(msg, outputFile);
        outputFile << std::endl;
        outputFile << "  Received signal: ";
        CodecUtils::outputMat(output_DLD.received_signal, outputFile);
        outputFile << std::endl;
        // of course the list size are going to be the maximum list size
        outputFile << "  List sizes searched: ";
        CodecUtils::outputMat(output_DLD.list_ranks, outputFile);
        outputFile << std::endl;
      }
      else {
        // if there is an error, we save the correct message, the incorrectly decoded message, and the received signal and record the list index.
        DLD_error++;
        number_of_errors++;
        outputFile << "Decoding Error: " << std::endl;
        outputFile << "  Transmitted message: ";
        CodecUtils::outputMat(msg, outputFile);
        outputFile << std::endl;
        outputFile << "  Received signal: ";
        CodecUtils::outputMat(output_DLD.received_signal, outputFile);
        outputFile << std::endl;
        outputFile << "  Decoded message: ";
        CodecUtils::outputMat(output_DLD.message, outputFile);
        outputFile << std::endl;
        // of course the list size are going to be the maximum list size
        outputFile << "  List sizes searched: ";
        CodecUtils::outputMat(output_DLD.list_ranks, outputFile);
        outputFile << std::endl;
      }

      // update list ranks
      for (int i = 0; i < expected_list_ranks.size(); ++i) {
        expected_list_ranks[i] += output_DLD.list_ranks[i];
      }

      // update in the outputfile
      if (number_of_trials % 2000 == 0) {
        outputFile << "Trial: " << number_of_trials << std::endl;
      }
    }

    ////////  DLD Decoding Statistics /////////
    for (int j = 0; j < expected_list_ranks.size(); ++j) {
      expected_list_ranks[j] /= number_of_trials;
    }
    std::cout << "DLD correct decoding: " << DLD_correct << " , Percentage: "
    << (double)DLD_correct/number_of_trials << std::endl; 
    std::cout << "DLD Exceeded list:" << DLD_list_exceeded << ", Percentage: "
    << (double)DLD_list_exceeded/number_of_trials << std::endl;
    std::cout << "DLD wrong decoding: " << DLD_error << " , Percentage: " << (double)DLD_error/number_of_trials
    << std::endl;

    std::cout << "expected list ranks: [ " << expected_list_ranks[0] << ", "
    << expected_list_ranks[1] << " ]" << std::endl;
    
    // distance spectrum statistics
    std::cout << "List ranks spectrum data: " << std::endl;
    std::cout << "  List 0: ";
    CodecUtils::printMat(DLD_list_0_size); std::cout << std::endl;
    std::cout << "  List 1: ";
    CodecUtils::printMat(DLD_list_1_size); std::cout << std::endl;


    DLD_correct_vec.push_back((double)DLD_correct/number_of_trials);
    DLD_list_exceed_vec.push_back((double)DLD_list_exceeded/number_of_trials);
    DLD_error_vec.push_back((double)DLD_error/number_of_trials);
  }


  ///////////  End of simulation ////////////
  std::cout << "DLD Decoder CORRECT decoding percentage: ";
  CodecUtils::print(DLD_correct_vec); std::cout << std::endl;
  std::cout << "DLD Decoder LIST EXCEEDED decoding percentage: ";
  CodecUtils::print(DLD_list_exceed_vec); std::cout << std::endl;
  std::cout << "DLD Decoder WRONG decoding percentage: ";
  CodecUtils::print(DLD_error_vec); std::cout << std::endl;
  outputFile.close();
}