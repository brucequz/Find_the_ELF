/**
 * @file main.cpp
 * @author Bruce Qu (brucequ@ucla.edu)
 * @brief 
 * @version 0.1
 * @date 2023-12-27
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <map>
#include <queue>
#include <random>
#include <sstream>
#include <chrono>

#include "../include/dualListDecoder.h"
#include "../include/viterbiCodec.h"
#include "../include/stopWatch.h"
// #include "mat.h"

static std::random_device rd{};
static std::mt19937 noise_gen(82);  // 82
int seed = 33; // 47
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

namespace {

struct mixed_info {
  std::vector<std::chrono::milliseconds> stepDurations;
  double DLD_decoded_portion;
  std::vector<int> DLD_expected_list_sizes;
};
}

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
void dualListExperiment(double snr_dB, int maximum_list_size,
                        int max_errors);
std::chrono::milliseconds softViterbiExperiment(double snr_dB, int max_errors);
mixed_info mixedDualListExperiment(double SNR_dB, int maximum_list_size, int max_errors);
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
  std::vector<double> SNR_dB;  // SNR is required for noise computations
  double offset = 10 * log10((double)2 * 64 /
                             (double)134);  // real rate of this code is 32/512
  for (int i = 0; i < EbN0.size(); i++) {
    SNR_dB.push_back(EbN0[i]);
  }
  std::cout << std::endl;
  std::cout << "Simulation begin: Running for snr points " << SNR_dB.front() << " -- " << SNR_dB.back() << std::endl;

  // 2. number of maximum errors to accumulate
  int max_errors = 50;
  int list_size = 1;

  // 3. DLD decoding setup
  // int DLD_maximum_list_size = 500;
  // std::vector<int> max_list_size_vector = {2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 1000};
  // std::vector<int> max_list_size_vector = {1};
  std::vector<int> max_list_size_vector = {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,500,750,1000};
  
  // step time recorder
  // mixed
  std::vector<std::chrono::milliseconds> mixed_ssv_time;
  std::vector<std::chrono::milliseconds> reconstruct_trellis_time;
  std::vector<std::chrono::milliseconds> mixed_insertion_time;
  std::vector<std::vector<int>> mixed_expected_list_ranks(2);
  std::vector<double> DLD_decoded_portion;
  // SSV
  std::vector<std::chrono::milliseconds> standard_soft_viterbi_time;
  for (double SNR : SNR_dB) {
    for (int DLD_maximum_list_size: max_list_size_vector) {
      std::cout << "Running experiments for max list size  = " << DLD_maximum_list_size << std::endl;
      mixed_info result = mixedDualListExperiment(SNR, DLD_maximum_list_size, max_errors);
      // record the results
      mixed_ssv_time.push_back(result.stepDurations[0]);
      reconstruct_trellis_time.push_back(result.stepDurations[1]);
      mixed_insertion_time.push_back(result.stepDurations[2]);
      mixed_expected_list_ranks[0].push_back(result.DLD_expected_list_sizes[0]);
      mixed_expected_list_ranks[1].push_back(result.DLD_expected_list_sizes[1]);
      DLD_decoded_portion.push_back(result.DLD_decoded_portion);
      // 4. soft Viterbi Experiment Setup
      
      std::chrono::milliseconds softTime = softViterbiExperiment(SNR, max_errors);
      standard_soft_viterbi_time.push_back(softTime);
    }
  }
  // time printing 
  std::cout << "printing mixed ssv time: [";
  for (std::chrono::milliseconds time : mixed_ssv_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "printing reconstruction trellis time: [";
  for (std::chrono::milliseconds time : reconstruct_trellis_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "printing mixed insertion time: [";
  for (std::chrono::milliseconds time : mixed_insertion_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "printing standard soft viterbi time: [";
  for (std::chrono::milliseconds time : standard_soft_viterbi_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "printing DLD 0 expected list ranks:";
  CodecUtils::print(mixed_expected_list_ranks[0]);

  std::cout << "printing DLD 1 expected list ranks:";
  CodecUtils::print(mixed_expected_list_ranks[1]);

  std::cout << "printing DLD decoded portion:";
  CodecUtils::print(DLD_decoded_portion);
  return 0;
}

std::chrono::milliseconds softViterbiExperiment(double snr_dB, int max_errors) {
  // output path
  // std::string outputFilePath = "../output/smaller_example/";
  // std::ofstream outputFile(outputFilePath + "STD_snr_" +
  //                          std::to_string(int(SNR_dB.front())) + "-" +
  //                          std::to_string(int(SNR_dB.back())) + ".txt");
  // if (!outputFile.is_open()) {
  //   std::cerr << "Failed to open the file for writing." << std::endl;
  // }

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
            << code.generator_poly[1] << " SNR = " << snr_dB << std::endl;
  
  // TODO: This is only for a single SNR_inplementation
  std::chrono::milliseconds softDurations(0);

  
  std::cout << "Now working on snr: " << snr_dB << "-------------------"
            << std::endl;

  std::vector<int> expected_list_ranks = {0, 0};

  int ML_decoding_error = 0;
  int Correct_decoding = 0;

  int number_of_trials = 0;
  int number_of_errors = 0;
  

  //while (number_of_errors < max_errors) {
  while (number_of_trials < 10000) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
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
    
    MessageInformation output_ML = codec.softViterbiDecoding(received_signal, softDurations);

    assert(output_ML.message.size() == msg.size());
    if (CodecUtils::areVectorsEqual(output_ML.message, msg)) {
      Correct_decoding++;
    } else {
      ML_decoding_error++;
      number_of_errors++;
    }

    // update in the outputfile
    // if (number_of_trials % 2000 == 0) {
    //   outputFile << "Trial: " << number_of_trials << std::endl;
    // }
  }

  std::cout << "For snr = " << snr_dB << ", " << number_of_trials
            << " were ran to accumulate " << number_of_errors << " errors."
            << std::endl;
  std::cout << "Time taken by soft viterbi trellis building and traceback: "
            << softDurations.count() << " milliseconds" << std::endl;
  std::cout << "std viterbi correct decoding: " << Correct_decoding
            << " , Percentage: "
            << (double)Correct_decoding / number_of_trials << std::endl;
  std::cout << "std viterbi wrong decoding: " << ML_decoding_error
            << " , Percentage: "
            << (double)ML_decoding_error / number_of_trials << std::endl;
  correct_decoding_snr.push_back((double)Correct_decoding / number_of_trials);
  ML_decoding_error_snr.push_back((double)ML_decoding_error /
                                  number_of_trials);

  ///////////  End of experiment ////////////
  std::cout << "Std viterbi CORRECT decoding percentage: ";
  CodecUtils::print(correct_decoding_snr);
  std::cout << std::endl;
  std::cout << "Std viterbi WRONG decoding percentage: ";
  CodecUtils::print(ML_decoding_error_snr);
  std::cout << std::endl;
  // outputFile.close();
  return softDurations;
}

mixed_info mixedDualListExperiment(double snr_dB, int maximum_list_size, int max_errors) {
  /*
  function description here
  TODO:

  */
  // output path
  // std::string outputFilePath = "../output/smaller_example/";
  // std::ostringstream oss;
  // oss << std::scientific << maximum_list_size;

  // std::string scientificString = oss.str();
  // std::ofstream outputFile(outputFilePath + "Mixed_" + scientificString + "_" +
  //                          std::to_string(int(SNR_dB.front())) + "-" +
  //                          std::to_string(int(SNR_dB.back())) + ".txt");
  // if (!outputFile.is_open()) {
  //   std::cerr << "Failed to open the file for writing." << std::endl;
  // }
  
  // return
  mixed_info result;

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

  // When using x^6+x^5+x^4+x+1 as generator polynomial, the octal representation is 163
  // Its decimal representation used for crc is 115

  // When using x^4+x^2+1 as generator polynomial, the octal representation is 25
  // Its decimal representation used for crc is 21

  // Therefore, the pair is either {21, 163} or {115, 25}
  code_1.k = 1;
  code_1.n = 1;
  code_1.v = 4;
  code_1.list_size = 1;
  code_1.crc_dec = 115;  // decimal
  code_1.crc_length = 7;  // highest degree + 1
  code_1.generator_poly = {25};  // octal

  CodeInformation code_2;
  // 56721
  // x^14+x^12+x^11+x^10+x^8+x^7+x^6+x^4+1 =
  //                         CRC: (x^2+x+1)
  //                         generator: x^12+x^11+x^10+x^9+x^8+x^5+x^3+x+1

  // 2473
  // x^10 + x^8 + x^5 + x^4 + x^3 + x + 1
  //                         CRC: (x^3+x+1)
  //                         generator: x^7+x^4+1

  // When using x^7+x^4+1 as generator polynomial, the octal representation is 221
  // Its decimal representation used for crc is 145

  // When using x^3+x+1 as generator polynomial, the octal representation is 13
  // Its decimal representation used for crc is 11

  // Therefore, the pair is either {11, 221} or {145, 13}
  code_2.k = 1;
  code_2.n = 1;
  code_2.v = 3;
  code_2.list_size = 1;
  code_2.crc_dec = 145;  // decimal
  code_2.crc_length = 8;  // highest degree + 1
  code_2.generator_poly = {13};  // octal 

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  std::vector<CodeInformation> dld_codes = {code_1, code_2};
  DualListDecoder DLD(dld_codes, maximum_list_size);

  // Simulation begin
  srand(seed);
  noise_gen.seed(82);
  int num_bits = 64;

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_list_exceed_vec;
  std::vector<double> DLD_error_vec;

  std::cout << "Running DLD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << snr_dB << std::endl;

  std::vector<std::chrono::milliseconds> stepTimeDurations(3, std::chrono::milliseconds(0));
  std::vector<int> expected_list_ranks = {1, 1};
  // vector to keep track of list sizes for both decoders
  std::vector<int> DLD_list_0_size;
  std::vector<int> DLD_list_1_size;

  int number_of_trials = 0;
  int number_of_errors = 0;

  std::cout << "Now working on snr: " << snr_dB << "-------------------"
            << std::endl;
  
  int SSV_correct = 0;
  int SSV_error = 0;
  int DLD_correct = 0;
  int DLD_list_exceeded = 0;
  int DLD_error = 0;

  std::vector<double> received_to_decoded_dist;
  std::vector<double> received_to_correct_dist;
  
  // measure the time of each operation
  // step1: SSV_list_1: time taken to add-compare-select and perform the initial trace-back
  // step2: C_trace_1: time taken for additional traaceback opearations required by SLVD
  // step3: C_insert_1: time taken of inserting new elements to maintain an ordered list of path metric differences

  //std::chrono::milliseconds softDurations(0);

  //while (number_of_errors < max_errors) {
  while (number_of_trials < 10000) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
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
    DLDInfo output_DLD = DLD.AdaptiveDecode(received_signal, stepTimeDurations);

    // if DLD declares List size exceeded
    // then resort to soft viterbi decoding
    if (output_DLD.message == std::vector<int>(64,-1)) {
      //std::cout << "List size EXCEEDED! Running SSV Now!" << std::endl;
      MessageInformation output_SSV = codec.softViterbiDecoding(received_signal, stepTimeDurations[1]);
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
    // if (number_of_trials % 2000 == 0) {
    //   outputFile << "Trial: " << number_of_trials << std::endl;
    // }
  }



  ////////  DLD Decoding Statistics /////////
  // record the time taken by each step
  for (size_t i = 0; i < stepTimeDurations.size(); ++i) {
    std::cout << "Step " << i + 1 << " Runtime: " << stepTimeDurations[i].count() << " milliseconds" << std::endl;
  }
  
  std::cout << "Out of " << number_of_trials << " experiments, " 
            << DLD_list_0_size.size() << " experiments were recorded to have DLD list sizes." << std::endl;
  

  // convert double to number with 2 digits precision
  std::ostringstream stream;
  stream << std::fixed << std::setprecision(2) << snr_dB;
  std::string snr_string = stream.str();
  
  std::string DLDList0_location = "../output/smaller_example/" + std::to_string(int(maximum_list_size)) + "_" + snr_string + "_list0.txt";
  std::ofstream outputDLDList0(DLDList0_location);
  if (outputDLDList0.is_open()) {
    for (const auto& element: DLD_list_0_size) {
      outputDLDList0 << element << "\n";
    }
    outputDLDList0.close();
    std::cout << "DLD List 0 successfully written!" << std::endl;
  } else {
    std::cerr << "Error opening the file!" << std::endl;
  }

  std::string DLDList1_location = "../output/smaller_example/" +std::to_string(int(maximum_list_size)) + "_" + snr_string + "_list1.txt";
  std::ofstream outputDLDList1(DLDList1_location);
  if (outputDLDList1.is_open()) {
    for (const auto& element: DLD_list_1_size) {
      outputDLDList1 << element << "\n";
    }
    outputDLDList1.close();
    std::cout << "DLD List 1 successfully written!" << std::endl;
  } else {
    std::cerr << "Error opening the file!" << std::endl;
  }
  
  std::cout << "Debug: DLD list 0 size: " << DLD_list_0_size.size() << std::endl; 
  for (int j = 0; j < expected_list_ranks.size(); ++j) {
    expected_list_ranks[j] /= DLD_list_0_size.size();
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
  // outputFile.close();
  

  // package result
  result.stepDurations = stepTimeDurations;
  result.DLD_decoded_portion = (double)DLD_list_0_size.size()/number_of_trials;
  result.DLD_expected_list_sizes = expected_list_ranks;
  return result;
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

  std::vector<std::chrono::milliseconds> timeDurations(3, std::chrono::milliseconds(0));

  DLDInfo output_DLD = DLD.AdaptiveDecode(received_signal, timeDurations);

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

    std::vector<std::chrono::milliseconds> timeDurations(3, std::chrono::milliseconds(0));
    std::chrono::milliseconds softDurations(0);

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
      DLDInfo output_DLD = DLD.AdaptiveDecode(received_signal, timeDurations);

      // if DLD declares List size exceeded
      // then resort to soft viterbi decoding
      MessageInformation output_SSV = codec.softViterbiDecoding(received_signal, softDurations);
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
            codec.softViterbiDecoding(received_signal, softDurations);
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
            codec.softViterbiDecoding(received_signal, softDurations);
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