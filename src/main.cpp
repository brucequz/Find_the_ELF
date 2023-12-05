#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <chrono>

#include "../include/dualListDecoder.h"
#include "../include/viterbiCodec.h"
// #include "mat.h"

static std::random_device rd{};
static std::mt19937 noise_gen(82);  // 82
int seed = 47;
static std::mt19937 msg_gen(seed);

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
// void writeToMat(const std::vector<std::vector<T>>& data, const char*
// filePath,
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
Functions that end with Experiment is a stand-alone experiment designed to run
by itself Functions taht end with Decoding is a decoding function that takes in
a single received message and
*/
void dualListExperiment(std::vector<double> SNR_dB, int maximum_list_size,
                        int max_errors);
void softViterbiExperiment(std::vector<double> SNR_dB, int max_errors);
void mixedDualListExperiment(std::vector<double> SNR_dB, int maximum_list_size, int max_errors);
bool dualListDecode(std::vector<double> received_signal,
                    std::vector<int> correct_message, int maximum_list_size);

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
  // int DLD_maximum_list_size = 500;
  std::vector<int> max_list_size_vector = {100, 200, 500, 1000};
  for (int DLD_maximum_list_size: max_list_size_vector) {
    std::cout << std::endl;
    std::cout << "Running experiments for max list size  = " << DLD_maximum_list_size << std::endl;
    auto DLD_start = std::chrono::high_resolution_clock::now();
    mixedDualListExperiment(SNR_dB, DLD_maximum_list_size, max_errors);
    auto DLD_end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(DLD_end - DLD_start);
    std::cout << "Time taken by function: " << duration.count() << " microseconds" << std::endl;

    // 4. soft Viterbi Experiment Setup
    auto SSV_start = std::chrono::high_resolution_clock::now();
    softViterbiExperiment(SNR_dB, max_errors);
    auto SSV_end = std::chrono::high_resolution_clock::now();
    auto SSV_duration = std::chrono::duration_cast<std::chrono::microseconds>(SSV_end - SSV_start);
    std::cout << "Time taken by function: " << SSV_duration.count() << " microseconds" << std::endl;
    std::cout << std::endl;
  }
  return 0;
}

void softViterbiExperiment(std::vector<double> SNR_dB, int max_errors) {
  // output path
  std::string outputFilePath = "../output/smaller_example/";
  std::ofstream outputFile(outputFilePath + "STD_snr_" +
                           std::to_string(int(SNR_dB.front())) + "-" +
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


  srand(seed);
  noise_gen.seed(82);
  int num_bits = 64;
  

  std::vector<double> correct_decoding_snr;
  std::vector<double> ML_decoding_error_snr;

  std::cout << "Running STD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << SNR_dB[0] << " - "
            << SNR_dB.back() << std::endl;

  for (double snr_dB : SNR_dB) {
    std::cout << "Now working on snr: " << snr_dB << "-------------------"
              << std::endl;

    std::vector<int> expected_list_ranks = {0, 0};

    int ML_decoding_error = 0;
    int Correct_decoding = 0;

    int number_of_trials = 0;
    int number_of_errors = 0;

    //while (number_of_errors < max_errors) {
    while (number_of_trials < 10) {
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
      // outputFile << " Printing generated message: " << std::endl;
      // CodecUtils::outputMat(msg, outputFile);
      // outputFile << std::endl;

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
  std::cout << "Std viterbi CORRECT decoding percentage: ";
  CodecUtils::print(correct_decoding_snr);
  std::cout << std::endl;
  std::cout << "Std viterbi WRONG decoding percentage: ";
  CodecUtils::print(ML_decoding_error_snr);
  std::cout << std::endl;
  outputFile.close();
}

void mixedDualListExperiment(std::vector<double> SNR_dB, int maximum_list_size, int max_errors) {
  /*
  function description here
  TODO:

  */
  // output path
  std::string outputFilePath = "../output/smaller_example/";
  std::ostringstream oss;
  oss << std::scientific << maximum_list_size;

  std::string scientificString = oss.str();
  std::ofstream outputFile(outputFilePath + "Mixed_" + scientificString + "_" +
                           std::to_string(int(SNR_dB.front())) + "-" +
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
  srand(seed);
  int num_bits = 64;

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_list_exceed_vec;
  std::vector<double> DLD_error_vec;

  std::cout << "Running DLD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << SNR_dB[0] << " - "
            << SNR_dB.back() << std::endl;

  for (double snr_dB : SNR_dB) {
    std::cout << "Now working on snr: " << snr_dB << "-------------------"
              << std::endl;

    std::vector<int> expected_list_ranks = {0, 0};
    
    int SSV_correct = 0;
    int SSV_error = 0;
    int DLD_correct = 0;
    int DLD_list_exceeded = 0;
    int DLD_error = 0;

    // vector to keep track of list sizes for both decoders
    std::vector<int> DLD_list_0_size;
    std::vector<int> DLD_list_1_size;

    int number_of_trials = 0;
    int number_of_errors = 0;

    std::vector<double> received_to_decoded_dist;
    std::vector<double> received_to_correct_dist;

    //while (number_of_errors < max_errors) {
    while (number_of_trials < 10) {
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
      // outputFile << " Printing generated message: " << std::endl;
      // CodecUtils::outputMat(msg, outputFile);
      // outputFile << std::endl;

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

      // if DLD declares List size exceeded
      // then resort to soft viterbi decoding
      if (output_DLD.message == std::vector<int>(64,-1)) {
        std::cout << "List size EXCEEDED! Running SSV Now!" << std::endl;
        MessageInformation output_SSV = codec.softViterbiDecoding(received_signal);
        if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
          SSV_correct++;
        } else {
          SSV_error++;
        }
        continue;
      }
      // we save the list ranks
      DLD_list_0_size.push_back(output_DLD.list_ranks[0]);
      DLD_list_1_size.push_back(output_DLD.list_ranks[1]);

      if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
        // CASE 1
        // correct decoding
        DLD_correct++;
      } else if (output_DLD.message == std::vector<int>(64, -1)) {
        // CASE 2
        // list size exceeded
        // we save the message, received_signal and record the maximum list size
        std::cerr << "THIS IS INCORRECT, UNANTICIPATED BEHAVIOR!!!" << std::endl;
      } else {
        // CASE 3
        // if there is an error, we save the correct message, the incorrectly
        // decoded message, and the received signal and record the list index.
        DLD_error++;
        number_of_errors++;
        std::cerr << "DLD decode find something, but it's incorrect" << std::endl;
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
    // record the distances
    
    std::cout << "Out of " << number_of_trials << " experiments, " 
              << DLD_list_0_size.size() << " experiments were recorded to have DLD list sizes." << std::endl;

    for (int j = 0; j < expected_list_ranks.size(); ++j) {
      expected_list_ranks[j] /= number_of_trials;
    }
    std::cout << "DLD correct decoding: " << DLD_correct
              << " , Percentage: " << (double)DLD_correct / number_of_trials
              << std::endl;
    std::cout << "DLD Exceeded list:" << DLD_list_exceeded << ", Percentage: "
              << (double)DLD_list_exceeded / number_of_trials << std::endl;
    std::cout << "DLD wrong decoding: " << DLD_error
              << " , Percentage: " << (double)DLD_error / number_of_trials
              << std::endl;
    std::cout << "SSV correct decoding: " << SSV_correct
              << " , Percentage: " << (double)SSV_correct / number_of_trials
              << std::endl;
    std::cout << "SSV wrong decoding: " << SSV_error
              << " , Percentage: " << (double)SSV_error / number_of_trials
              << std::endl;

    std::cout << "Expected list ranks: [ " << expected_list_ranks[0] << ", "
              << expected_list_ranks[1] << " ]" << std::endl;

    DLD_correct_vec.push_back((double)(DLD_correct+SSV_correct) / number_of_trials);
    DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
    DLD_error_vec.push_back((double)DLD_error / number_of_trials);
  }

  ///////////  End of simulation ////////////
  std::cout << "Mixed Decoder CORRECT decoding percentage: ";
  CodecUtils::print(DLD_correct_vec);
  std::cout << std::endl;
  std::cout << "Mixed Decoder LIST EXCEEDED decoding percentage: ";
  CodecUtils::print(DLD_list_exceed_vec);
  std::cout << std::endl;
  std::cout << "Mixed Decoder WRONG decoding percentage: ";
  CodecUtils::print(DLD_error_vec);
  std::cout << std::endl;
  outputFile.close();
}

bool dualListDecode(std::vector<double> received_signal,
                    std::vector<int> correct_message, int maximum_list_size) {
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

  int DLD_list_0_size;
  int DLD_list_1_size;

  DLDInfo output_DLD = DLD.adaptiveDecode(received_signal);

  DLD_list_0_size = output_DLD.list_ranks[0];
  DLD_list_1_size = output_DLD.list_ranks[1];

  std::cout << maximum_list_size << " DLD decoding combined path metric: " << output_DLD.combined_metric << std::endl;

  if (CodecUtils::areVectorsEqual(output_DLD.message, correct_message)) {
    // CASE 1
    // correct decoding
    std::cout << "  Found correct decoding at ";
    CodecUtils::printMat(output_DLD.list_ranks); std::cout << std::endl;
    CodecUtils::printMat(output_DLD.message); std::cout << std::endl;
    CodecUtils::printMat(correct_message); std::cout << std::endl;
    std::cout << std::endl;
    return true;
  } else if (output_DLD.message == std::vector<int>(64, -1)) {
    // CASE 2
    // list size exceeded
    // we save the message, received_signal and record the maximum list size
    std::cout << "LSE Error: " << std::endl;
    return false;
  } else {
    // CASE 3
    // if there is an error, we save the correct message, the incorrectly
    // decoded message, and the received signal and record the list index.
    std::cout << "Decoding Error: " << std::endl;
    return false;
  }
}

void dualListExperiment(std::vector<double> SNR_dB, int maximum_list_size,
                        int max_errors) {
  /*
  function description here
  TODO:

  */
  // output path
  std::string outputFilePath = "../output/smaller_example/";
  std::ostringstream oss;
  oss << std::scientific << maximum_list_size;

  std::string scientificString = oss.str();
  std::ofstream outputFile(outputFilePath + "DLD_" + scientificString + "_" +
                           std::to_string(int(SNR_dB.front())) + "-" +
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

  std::cout << "Running DLD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << SNR_dB[0] << " - "
            << SNR_dB.back() << std::endl;

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

    std::vector<double> received_to_decoded_dist;
    std::vector<double> received_to_correct_dist;

    //while (number_of_errors < max_errors) {
    while (number_of_trials < 10000) {
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

      // msg = {0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 0, 0,
      //        1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0, 1, 0,
      //        0, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 1, 0, 1, 0, 0, 0};

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
      // received_signal = {
      //     0.933018,  0.524545,  1.96875,   1.82899,     0.297545,  1.28138,
      //     -1.41533,  -0.766253, -0.583526, 0.000501252, 1.44843,   1.34277,
      //     0.434987,  0.614327,  -1.9752,   -0.172404,   -1.34518,  0.982204,
      //     -2.1783,   -1.99824,  -0.247085, -1.34994,    0.0391971, -1.59384,
      //     1.3158,    -0.878993, 1.2844,    1.03464,     0.328343,  -0.735207,
      //     1.23869,   -1.72009,  -1.29932,  -0.47016,    -1.79555,  0.910431,
      //     1.74748,   -1.74393,  -1.33471,  -1.13221,    0.448533,  -0.826694,
      //     0.679492,  -1.50642,  -0.13838,  0.965305,    0.377865,  -0.287025,
      //     0.488627,  0.563657,  -1.3423,   -0.576363,   1.67196,   1.05073,
      //     1.68197,   0.597487,  1.32816,   -1.21886,    -0.892688, -0.0229954,
      //     1.19326,   0.666117,  -2.04728,  -1.6684,     0.24426,   0.395314,
      //     -0.444505, -1.71933,  1.25274,   0.906817,    -1.5019,   -1.22162,
      //     1.46852,   -1.42066,  -0.461029, -0.943509,   -0.978827, 1.12364,
      //     -0.329824, 1.62396,   -1.44469,  0.322167,    1.27247,   0.683425,
      //     1.00047,   -0.80628,  -1.58197,  -0.221272,   1.25508,   1.12334,
      //     0.796071,  -0.65488,  -1.01318,  0.1393,      -1.5108,   -0.66198,
      //     -1.28406,  -1.28309,  2.1364,    -0.832802,   -0.832442, 0.581645,
      //     0.67589,   1.21626,   1.5545,    -0.145063,   0.722937,  0.103384,
      //     -0.688162, 1.79378,   -1.77359,  0.323369,    -0.726902, 0.458736,
      //     -1.8405,   0.308128,  1.14545,   -1.18301,    1.51812,   0.93853,
      //     -0.955463, 0.0436563, -1.26439,  1.15478,     -0.693205, 2.15561,
      //     0.454222,  -0.472711, -1.93238,  1.93998,     1.46456,   0.798949,
      //     1.61938,   -0.527419, -1.30957,  -1.34847,    -0.289711, 0.215369,
      //     -1.46554,  -0.145958, -0.38653,  -0.531511,   1.73939,   0.353206,
      //     0.0228116, 0.22973,   2.44116,   1.02029};

      // DLD DECODING
      DLDInfo output_DLD = DLD.adaptiveDecode(received_signal);

      // if DLD declares List size exceeded
      // then resort to soft viterbi decoding
      MessageInformation output_SSV = codec.softViterbiDecoding(received_signal);
      // we save the list ranks
      DLD_list_0_size.push_back(output_DLD.list_ranks[0]);
      DLD_list_1_size.push_back(output_DLD.list_ranks[1]);

      if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
        // CASE 1
        // correct decoding
        // record the distance between the decoded and the correct
        received_to_correct_dist.push_back(output_DLD.combined_metric);
        received_to_decoded_dist.push_back(output_DLD.combined_metric);
        DLD_correct++;
      } else if (output_DLD.message == std::vector<int>(64, -1)) {
        // CASE 2
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

        MessageInformation ssv_output =
            codec.softViterbiDecoding(received_signal);
        if (ssv_output.message == msg) {
          
          outputFile << " SSV decoded successfully!!!" << std::endl;
          outputFile << " Now running DLD with a larger list size" << std::endl;
          bool decoding_5e6 = dualListDecode(received_signal, msg, 5e6);
          outputFile << " 5e6 list size Decoding result: "
                     <<  decoding_5e6 << std::endl;
          // outputFile << " 1e7 list size Decoding result: "
                    //  << dualListDecode(received_signal, msg, 1e7) << std::endl;
          // record the distance
          received_to_correct_dist.push_back(ssv_output.path_metric);
          if (decoding_5e6 == true) {
            received_to_decoded_dist.push_back(ssv_output.path_metric);
            std::cout << "SSV decoding path metric " << ssv_output.path_metric << std::endl; 
          } else {
            received_to_decoded_dist.push_back(-1.0);
          }
        } else {
          outputFile << " SSV decoding insuccessful" << std::endl;
          received_to_correct_dist.push_back(-1.0);
          received_to_decoded_dist.push_back(-1.0);
        }
        number_of_errors = max_errors + 1;
        break;
      } else {
        // CASE 3
        // if there is an error, we save the correct message, the incorrectly
        // decoded message, and the received signal and record the list index.
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
        outputFile << "  List sizes searched: ";
        CodecUtils::outputMat(output_DLD.list_ranks, outputFile);
        outputFile << std::endl;

        MessageInformation ssv_output =
            codec.softViterbiDecoding(received_signal);
        if (ssv_output.message == msg) {
          outputFile << " SSV decoded successfully!!!" << std::endl;
          received_to_correct_dist.push_back(ssv_output.path_metric);
          received_to_decoded_dist.push_back(output_DLD.combined_metric);
        } else {
          outputFile << " SSV decoding insuccessful" << std::endl;
          received_to_correct_dist.push_back(-1.0);
          received_to_decoded_dist.push_back(output_DLD.combined_metric);
        }
        number_of_errors = max_errors + 1;
        break;
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
    // record the distances
    assert(received_to_correct_dist.size() == 1);
    outputFile << "Printing distances to the correct codeword " << std::endl;
    dualdecoderutils::outputMat(received_to_correct_dist, outputFile);
    outputFile << std::endl;
    
    assert(received_to_decoded_dist.size() == 1);
    outputFile << "Printing distances to the decoded codeword " << std::endl;
    dualdecoderutils::outputMat(received_to_decoded_dist, outputFile);
    outputFile << std::endl;

    for (int j = 0; j < expected_list_ranks.size(); ++j) {
      expected_list_ranks[j] /= number_of_trials;
    }
    std::cout << "DLD correct decoding: " << DLD_correct
              << " , Percentage: " << (double)DLD_correct / number_of_trials
              << std::endl;
    std::cout << "DLD Exceeded list:" << DLD_list_exceeded << ", Percentage: "
              << (double)DLD_list_exceeded / number_of_trials << std::endl;
    std::cout << "DLD wrong decoding: " << DLD_error
              << " , Percentage: " << (double)DLD_error / number_of_trials
              << std::endl;

    std::cout << "expected list ranks: [ " << expected_list_ranks[0] << ", "
              << expected_list_ranks[1] << " ]" << std::endl;

    // distance spectrum statistics
    std::cout << "List ranks spectrum data: " << std::endl;
    std::cout << "  List 0: ";
    CodecUtils::printMat(DLD_list_0_size);
    std::cout << std::endl;
    std::cout << "  List 1: ";
    CodecUtils::printMat(DLD_list_1_size);
    std::cout << std::endl;

    DLD_correct_vec.push_back((double)DLD_correct / number_of_trials);
    DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
    DLD_error_vec.push_back((double)DLD_error / number_of_trials);
  }

  ///////////  End of simulation ////////////
  std::cout << "DLD Decoder CORRECT decoding percentage: ";
  CodecUtils::print(DLD_correct_vec);
  std::cout << std::endl;
  std::cout << "DLD Decoder LIST EXCEEDED decoding percentage: ";
  CodecUtils::print(DLD_list_exceed_vec);
  std::cout << std::endl;
  std::cout << "DLD Decoder WRONG decoding percentage: ";
  CodecUtils::print(DLD_error_vec);
  std::cout << std::endl;
  outputFile.close();
}