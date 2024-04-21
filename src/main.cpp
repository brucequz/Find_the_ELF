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
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <random>
#include <sstream>

#include "../include/dualListDecoder.h"
#include "../include/stopWatch.h"
#include "../include/viterbiCodec.h"
#include "../include/CONSTANTS.h"
// #include "mat.h"

static std::random_device rd{};
static std::mt19937 noise_gen(82);  // 82
int seed = 33;                      // 47
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

template <typename T>
void print(const std::vector<T>& vec) {
  for (const T& element : vec) {
    std::cout << element << " ";
  }
  std::cout << " with size: " << vec.size();
  std::cout << std::endl;
}

template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      std::cout << element << " ";
    }
    std::cout << ";" << std::endl;
  }
}
}  // namespace

/*
Functions that end with Experiment is a stand-alone experiment designed to run
by itself Functions taht end with Decoding is a decoding function that takes in
a single received message and
*/
void TimeAndComplexitySimulation(double SNR);
std::chrono::milliseconds softViterbiExperiment(double snr_dB, int max_errors);
mixed_info mixedDualListExperiment(double SNR_dB, int maximum_list_size,
                                   int max_errors);
bool dualListDecode(std::vector<double> received_signal,
                    std::vector<int> correct_message, int maximum_list_size);
void dualListExperiment_rate_1_3(double snr_dB, int max_list_size,
                                 int max_errors);
mixed_info mixedDualListExperiment_rate_1_3(double snr_dB,
                                            int maximum_list_size,
                                            int max_errors);
std::chrono::milliseconds softViterbiExperiment_rate_1_3(double snr_dB,
                                                         int max_errors);
void TimeAndComplexitySimulation_rate_1_3(double SNR);

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
  double offset;
  if (PUNCTURE) {
    offset = 10 * log10((double)2 * K /
                        (double)(PUNC_RATE *
                                 (K + V)));  // real rate of this code is 32/512
  } else {
    offset =
        10 *
        log10((double)2 * K /
              (double)(RDEN * (K + V)));  // real rate of this code is 32/512
  }
  std::cout << "offset : " << offset << std::endl;
  for (int i = 0; i < EbN0.size(); i++) {
    SNR_dB.push_back(EbN0[i] + offset);
  }

  ////////////////////// Dual List Decoding Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   dualListExperiment(snr_dB, 1e5, 1000);
  // }

  ////////////////////// Time and Complexity Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   TimeAndComplexitySimulation(snr_dB);
  // }

  ////////////////////// DSU/Soft Viterbi Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   softViterbiExperiment(snr_dB, 5);
  // }

  ////////////////////// DSU/Mixed Decoder Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   mixedDualListExperiment(snr_dB, 300, 200);
  // }

  ////////////////////// RATE 1/3 Dual List Decoding Experiment
  /////////////////////////
  for (double snr_dB : SNR_dB) {
    dualListExperiment_rate_1_3(snr_dB, MAX_LIST_SIZE, 20);
  }

  ////////////////////// RATE 1/3 DSU/Mixed Decoder Experiment
  /////////////////////////
  // for (double snr_dB : SNR_dB) {
  //   mixedDualListExperiment_rate_1_3(snr_dB, MAX_LIST_SIZE, 50);
  // }

  ////////////////////// RATE 1/3 DSU/Soft Viterbi Experiment
  /////////////////////////
  // for (double snr_dB : SNR_dB) {
  //   softViterbiExperiment_rate_1_3(snr_dB, 5);
  // }

  ////////////////////// Time and Complexity Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   TimeAndComplexitySimulation_rate_1_3(snr_dB);
  // }

  ////////////////////// Find D_free //////////////////////
  // FindDFree_ADFree();
  return 0;
}

void TimeAndComplexitySimulation(double SNR) {
  std::cout << std::endl;
  std::cout << "Simulation begin: Running for snr points " << SNR << std::endl;

  // 2. number of maximum errors to accumulate
  int max_errors = 50;
  int list_size = 1;

  // 3. DLD decoding setup
  // int DLD_maximum_list_size = 500;
  std::vector<int> max_list_size_vector = {2,   5,   10,  20,  50,  100,
                                           200, 300, 400, 500, 1000};
  // std::vector<int> max_list_size_vector = {1,2,3,4,5,6,7,8,9,10,12,14,16,18};
  // std::vector<int> max_list_size_vector = {70};
  // std::vector<int> max_list_size_vector =
  // {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,500,750,1000};
  // step time recorder
  // mixed
  std::vector<std::chrono::milliseconds> mixed_ssv_time;
  std::vector<std::chrono::milliseconds> reconstruct_trellis_time;
  std::vector<std::chrono::milliseconds> mixed_insertion_time;
  std::vector<std::vector<int>> mixed_expected_list_ranks(2);
  std::vector<double> DLD_decoded_portion;
  // SSV
  std::vector<std::chrono::milliseconds> standard_soft_viterbi_time;
  for (int DLD_maximum_list_size : max_list_size_vector) {
    std::cout << "Running experiments for max list size  = "
              << DLD_maximum_list_size << std::endl;
    mixed_info result =
        mixedDualListExperiment(SNR, DLD_maximum_list_size, max_errors);
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
}

void TimeAndComplexitySimulation_rate_1_3(double SNR) {
  std::cout << std::endl;
  std::cout << "Simulation begin: Running for snr points " << SNR << std::endl;

  // 2. number of maximum errors to accumulate
  int max_errors = 50;
  int list_size = 1;

  // 3. DLD decoding setup
  // int DLD_maximum_list_size = 500;
  // std::vector<int> max_list_size_vector = {2, 5, 10, 20, 50, 100, 200, 300,
  // 400, 500, 1000};
  // std::vector<int> max_list_size_vector =
  // {1100, 1200, 1300, 1400, 1500, 1800, 2000, 2200, 2300, 2500, 3000};
  std::vector<int> max_list_size_vector = {MAX_LIST_SIZE};
  // std::vector<int> max_list_size_vector =
  // {1,2,3,4,5,6,7,8,9,10,12,14,16,18,20,25,30,35,40,45,50,60,70,80,90,100,150,200,250,500,750,1000};
  // step time recorder
  // mixed
  std::vector<std::chrono::milliseconds> mixed_ssv_time;
  std::vector<std::chrono::milliseconds> reconstruct_trellis_time;
  std::vector<std::chrono::milliseconds> mixed_insertion_time;
  std::vector<std::vector<int>> mixed_expected_list_ranks(2);
  std::vector<double> DLD_decoded_portion;
  // SSV
  std::vector<std::chrono::milliseconds> standard_soft_viterbi_time;
  for (int DLD_maximum_list_size : max_list_size_vector) {
    std::cout << "Running experiments for max list size  = "
              << DLD_maximum_list_size << std::endl;
    mixed_info result = mixedDualListExperiment_rate_1_3(
        SNR, DLD_maximum_list_size, max_errors);
    // record the results
    mixed_ssv_time.push_back(result.stepDurations[0]);
    reconstruct_trellis_time.push_back(result.stepDurations[1]);
    mixed_insertion_time.push_back(result.stepDurations[2]);
    mixed_expected_list_ranks[0].push_back(result.DLD_expected_list_sizes[0]);
    mixed_expected_list_ranks[1].push_back(result.DLD_expected_list_sizes[1]);
    DLD_decoded_portion.push_back(result.DLD_decoded_portion);
    // 4. soft Viterbi Experiment Setup

    std::chrono::milliseconds softTime =
        softViterbiExperiment_rate_1_3(SNR, max_errors);
    standard_soft_viterbi_time.push_back(softTime);
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
  code.v = 9;
  code.list_size = 1;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {1157, 1753};  // octal
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

  // while (number_of_errors < max_errors) {
  while (number_of_trials < TRIALS) {
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

    MessageInformation output_ML =
        codec.softViterbiDecoding(received_signal, softDurations);

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
            << " , Percentage: " << (double)Correct_decoding / number_of_trials
            << std::endl;
  std::cout << "std viterbi wrong decoding: " << ML_decoding_error
            << " , Percentage: " << (double)ML_decoding_error / number_of_trials
            << std::endl;
  correct_decoding_snr.push_back((double)Correct_decoding / number_of_trials);
  ML_decoding_error_snr.push_back((double)ML_decoding_error / number_of_trials);

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
  code.generator_poly = {3345, 3613};  // octal
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
  code_1.generator_poly = {147};  // octal

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
  code_2.crc_dec = 13;
  code_2.crc_length = 4;
  code_2.generator_poly = {211};

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
    DLDInfo output_DLD = DLD.AdaptiveDecode_SimpleAlternate(received_signal, stepTimeDurations);

    // if DLD declares List size exceeded
    // then resort to soft viterbi decoding
    if (output_DLD.message == std::vector<int>(64,-1)) {
      //std::cout << "List size EXCEEDED! Running SSV Now!" << std::endl;
      MessageInformation output_SSV = codec.softViterbiDecoding(received_signal, stepTimeDurations[1]);
      if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
        SSV_correct++;
      } else {
        SSV_error++;
        number_of_errors++;
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


// Mixed New
// mixed_info mixedDualListExperiment(double snr_dB, int maximum_list_size,
//                                    int max_errors) {
//   /*
//   function description here
//   TODO:

//   */
//   // output path
//   // std::string outputFilePath = "../output/smaller_example/";
//   // std::ostringstream oss;
//   // oss << std::scientific << maximum_list_size;

//   // std::string scientificString = oss.str();
//   // std::ofstream outputFile(outputFilePath + "Mixed_" + scientificString + "_"
//   // +
//   //                          std::to_string(int(SNR_dB.front())) + "-" +
//   //                          std::to_string(int(SNR_dB.back())) + ".txt");
//   // if (!outputFile.is_open()) {
//   //   std::cerr << "Failed to open the file for writing." << std::endl;
//   // }

//   // return
//   mixed_info result;

//   /// Building Codec
//   CodeInformation code;
//   code.k = 1;
//   code.n = 2;
//   code.v = 10;
//   code.list_size = 1;
//   code.crc_dec = -1;
//   code.crc_length = -1;
//   code.generator_poly = {3345,
//                          3613};  // octal // might want to try {3345 and 3613}
//   ViterbiCodec codec(code);

//   CodeInformation code_1;
//   // 61713
//   // x^14+x^13+x^9+x^8+x^7+x^6+x^3+x+1 =
//   //                        CRC: (x^3+x^2+1)
//   //                        generator: (x^11+x^8+x^7+x^3+x^2+x+1)

//   // 3217: x^10 + x^9 + x^7 + x^3 + x^2 + x + 1
//   //                        CRC: (x^4+x^2+1)
//   //                        generator: x^6+x^5+x^4+x+1

//   // 3613: 11110001011 x^10 + x^9 + x^8 + x^7 + x^3 + x + 1
//   // One can divide it into: (x^6 + x^5 + x^2 + x + 1)*(x^4 + x^2 + 1)
//   // x^6 + x^5 + x^2 + x + 1 is (binary)1100111, (octal)147, (decimal)103
//   // x^4 + x^2 + 1 is (binary)10101, (octal)25, (decimal)21

//   code_1.k = 1;
//   code_1.n = 1;
//   code_1.v = 6;
//   code_1.list_size = 1;
//   code_1.crc_dec = 21;            // decimal
//   code_1.crc_length = 5;          // highest degree + 1
//   code_1.generator_poly = {147};  // octal

//   CodeInformation code_2;
//   // 56721
//   // x^14+x^12+x^11+x^10+x^8+x^7+x^6+x^4+1 =
//   //                         CRC: (x^2+x+1)
//   //                         generator: x^12+x^11+x^10+x^9+x^8+x^5+x^3+x+1

//   // 2473
//   // x^10 + x^8 + x^5 + x^4 + x^3 + x + 1
//   //                         CRC: (x^3+x+1)
//   //                         generator: x^7+x^4+1

//   // 3345: 11011100101 x^10 + x^9 + x^7 + x^6 + x^5 + x^2 + 1
//   // One can divide it into: (x^3 + x^2 + 1)*(x^7 + x^3 + 1)
//   // x^7 + x^3 + 1 is (binary)10001001, (octal)211, (decimal)137
//   // x^3 + x^2 + 1 is (binary)1101, (octal)15, (decimal)13

//   code_2.k = 1;
//   code_2.n = 1;
//   code_2.v = 7;
//   code_2.list_size = 1;
//   code_2.crc_dec = 13;            // decimal
//   code_2.crc_length = 4;          // highest degree + 1
//   code_2.generator_poly = {211};  // octal

//   std::vector<CodeInformation> dld_codes = {code_1, code_2};
//   DualListDecoder DLD(dld_codes, maximum_list_size);

//   // Simulation begin
//   srand(seed);
//   noise_gen.seed(82);
//   int num_bits = 64;

//   std::vector<double> DLD_correct_vec;
//   std::vector<double> DLD_list_exceed_vec;
//   std::vector<double> DLD_error_vec;

//   std::cout << "Running DLD tests for v = " << code.v
//             << " test, with generator poly: " << code.generator_poly[0] << ", "
//             << code.generator_poly[1] << " SNR = " << snr_dB << std::endl;

//   std::cout << "DLD order: g1 = " << code_1.generator_poly[0]
//             << ", g2 = " << code_2.generator_poly[0] << std::endl;

//   std::vector<std::chrono::milliseconds> stepTimeDurations(
//       3, std::chrono::milliseconds(0));
//   std::vector<int> expected_list_ranks = {1, 1};
//   // vector to keep track of list sizes for both decoders
//   std::vector<int> DLD_list_0_size;
//   std::vector<int> DLD_list_1_size;

//   int number_of_trials = 0;
//   int number_of_errors = 0;

//   std::cout << "Now working on snr: " << snr_dB << "-------------------"
//             << std::endl;

//   int SSV_correct = 0;
//   int SSV_error = 0;
//   int DLD_correct = 0;
//   int DLD_list_exceeded = 0;
//   int DLD_error = 0;

//   std::vector<double> received_to_decoded_dist;
//   std::vector<double> received_to_correct_dist;

//   // measure the time of each operation
//   // step1: SSV_list_1: time taken to add-compare-select and perform the initial
//   // trace-back step2: C_trace_1: time taken for additional traaceback
//   // opearations required by SLVD step3: C_insert_1: time taken of inserting new
//   // elements to maintain an ordered list of path metric differences

//   // std::chrono::milliseconds softDurations(0);

//   // while (number_of_errors < max_errors) {
//   while (number_of_trials < TRIALS) {
//     number_of_trials++;

//     if (number_of_trials % 2000 == 0) {
//       std::cout << "Trial number: " << number_of_trials << std::endl;
//       std::cout << "Current number of errors: " << number_of_errors
//                 << std::endl;
//     }

//     std::vector<int> msg;
//     for (int i = 0; i < num_bits; ++i) {
//       int random_bit = rand() % 2;
//       msg.push_back(random_bit);
//     }
//     // outputFile << " Printing generated message: " << std::endl;
//     // CodecUtils::outputMat(msg, outputFile);
//     // outputFile << std::endl;

//     // coding
//     std::vector<int> encoded_msg = codec.encodeZTCC(msg);
//     // outputFile << " Printing coded message: " << std::endl;
//     // CodecUtils::outputMat(encoded_msg, outputFile);
//     // outputFile << std::endl;

//     assert(encoded_msg.size() == (msg.size() + code.v) * 2);

//     std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
//     // outputFile << "Printing modulated signal: " << std::endl;
//     // CodecUtils::outputMat(modulated_signal, outputFile);
//     // outputFile << std::endl;

//     std::vector<double> received_signal =
//         AWGN::addNoise(modulated_signal, snr_dB);
//     // outputFile << "Printing received signal: " << std::endl;
//     // CodecUtils::outputMat(received_signal, outputFile);
//     // outputFile << std::endl;

//     // TODO pucturing

//     // DLD DECODING
//     DLDInfo output_DLD =
//         DLD.adaptiveDecode(received_signal, stepTimeDurations);

//     // if DLD declares List size exceeded
//     // then resort to soft viterbi decoding
//     if (output_DLD.message == std::vector<int>(64, -1)) {
//       // std::cout << "List size EXCEEDED! Running SSV Now!" << std::endl;
//       MessageInformation output_SSV =
//           codec.softViterbiDecoding(received_signal, stepTimeDurations[1]);
//       if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
//         SSV_correct++;
//       } else {
//         SSV_error++;
//         number_of_errors++;
//         // std::cerr << "SSV decode find something, but it's incorrect" <<
//         // std::endl;
//       }
//       continue;
//     }
//     // we save the list ranks
//     DLD_list_0_size.push_back(output_DLD.list_ranks[0]);
//     DLD_list_1_size.push_back(output_DLD.list_ranks[1]);

//     if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
//       // CASE 1
//       // correct decoding
//       DLD_correct++;
//     } else if (output_DLD.message == std::vector<int>(64, -1)) {
//       // CASE 2
//       // list size exceeded
//       // we save the message, received_signal and record the maximum list size
//       std::cerr << "THIS IS INCORRECT, UNANTICIPATED BEHAVIOR!!!" << std::endl;
//     } else {
//       // CASE 3
//       // if there is an error, we save the correct message, the incorrectly
//       // decoded message, and the received signal and record the list index.
//       DLD_error++;
//       number_of_errors++;
//       std::cerr << "DLD decode find something, but it's incorrect" << std::endl;
//     }

//     // update list ranks
//     for (int i = 0; i < expected_list_ranks.size(); ++i) {
//       expected_list_ranks[i] += output_DLD.list_ranks[i];
//     }

//     // update in the outputfile
//     // if (number_of_trials % 2000 == 0) {
//     //   outputFile << "Trial: " << number_of_trials << std::endl;
//     // }
//   }

//   ////////  DLD Decoding Statistics /////////
//   // record the time taken by each step
//   for (size_t i = 0; i < stepTimeDurations.size(); ++i) {
//     std::cout << "Step " << i + 1
//               << " Runtime: " << stepTimeDurations[i].count() << " milliseconds"
//               << std::endl;
//   }

//   std::cout << "Out of " << number_of_trials << " experiments, "
//             << DLD_list_0_size.size()
//             << " experiments were recorded to have DLD list sizes."
//             << std::endl;

//   // convert double to number with 2 digits precision
//   std::ostringstream stream;
//   stream << std::fixed << std::setprecision(2) << snr_dB;
//   std::string snr_string = stream.str();

//   std::string DLDList0_location = "../output/smaller_example/" +
//                                   std::to_string(int(maximum_list_size)) + "_" +
//                                   snr_string + "_list0.txt";
//   std::ofstream outputDLDList0(DLDList0_location);
//   if (outputDLDList0.is_open()) {
//     for (const auto& element : DLD_list_0_size) {
//       outputDLDList0 << element << "\n";
//     }
//     outputDLDList0.close();
//     std::cout << "DLD List 0 successfully written!" << std::endl;
//   } else {
//     std::cerr << "Error opening the file!" << std::endl;
//   }

//   std::string DLDList1_location = "../output/smaller_example/" +
//                                   std::to_string(int(maximum_list_size)) + "_" +
//                                   snr_string + "_list1.txt";
//   std::ofstream outputDLDList1(DLDList1_location);
//   if (outputDLDList1.is_open()) {
//     for (const auto& element : DLD_list_1_size) {
//       outputDLDList1 << element << "\n";
//     }
//     outputDLDList1.close();
//     std::cout << "DLD List 1 successfully written!" << std::endl;
//   } else {
//     std::cerr << "Error opening the file!" << std::endl;
//   }

//   std::cout << "Debug: DLD list 0 size: " << DLD_list_0_size.size()
//             << std::endl;
//   for (int j = 0; j < expected_list_ranks.size(); ++j) {
//     expected_list_ranks[j] /= DLD_list_0_size.size();
//   }
//   std::cout << "DLD correct decoding: " << DLD_correct
//             << " , Percentage: " << (double)DLD_correct / number_of_trials
//             << std::endl;
//   std::cout << "DLD Exceeded list:" << DLD_list_exceeded
//             << ", Percentage: " << (double)DLD_list_exceeded / number_of_trials
//             << std::endl;
//   std::cout << "DLD wrong decoding: " << DLD_error
//             << " , Percentage: " << (double)DLD_error / number_of_trials
//             << std::endl;
//   std::cout << "SSV correct decoding: " << SSV_correct
//             << " , Percentage: " << (double)SSV_correct / number_of_trials
//             << std::endl;
//   std::cout << "SSV wrong decoding: " << SSV_error
//             << " , Percentage: " << (double)SSV_error / number_of_trials
//             << std::endl;

//   std::cout << "Expected list ranks: [ " << expected_list_ranks[0] << ", "
//             << expected_list_ranks[1] << " ]" << std::endl;

//   DLD_correct_vec.push_back((double)(DLD_correct + SSV_correct) /
//                             number_of_trials);
//   DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
//   DLD_error_vec.push_back((double)DLD_error + SSV_error / number_of_trials);

//   ///////////  End of simulation ////////////
//   std::cout << "Mixed Decoder CORRECT decoding percentage: ";
//   CodecUtils::print(DLD_correct_vec);
//   std::cout << std::endl;
//   std::cout << "Mixed Decoder LIST EXCEEDED decoding percentage: ";
//   CodecUtils::print(DLD_list_exceed_vec);
//   std::cout << std::endl;
//   std::cout << "Mixed Decoder WRONG decoding percentage: ";
//   CodecUtils::print(DLD_error_vec);
//   std::cout << std::endl;
//   // outputFile.close();

//   // package result
//   result.stepDurations = stepTimeDurations;
//   result.DLD_decoded_portion =
//       (double)DLD_list_0_size.size() / number_of_trials;
//   result.DLD_expected_list_sizes = expected_list_ranks;
//   return result;
// }

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
  //                         generator: x^7+x^4+1 (221 in octal)
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

  std::vector<std::chrono::milliseconds> timeDurations(
      3, std::chrono::milliseconds(0));

  DLDInfo output_DLD =
      DLD.AdaptiveDecode_SimpleAlternate(received_signal, timeDurations);

  DLD_list_0_size = output_DLD.list_ranks[0];
  DLD_list_1_size = output_DLD.list_ranks[1];

  std::cout << "Maximum list size: " << maximum_list_size
            << ", DLD decoding combined path metric: "
            << output_DLD.combined_metric << std::endl;

  if (CodecUtils::areVectorsEqual(output_DLD.message, correct_message)) {
    // CASE 1
    // correct decoding
    std::cout << "  Found correct decoding at ";
    CodecUtils::printMat(output_DLD.list_ranks);
    std::cout << std::endl;
    CodecUtils::printMat(output_DLD.message);
    std::cout << std::endl;
    CodecUtils::printMat(correct_message);
    std::cout << std::endl;
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

void dualListExperiment_rate_1_3(double snr_dB, int max_list_size,
                                 int max_errors) {
  /*
  function description here
  TODO:

  */
  // output path
  // std::string outputFilePath = "../output/list_size_exceeded/";
  // std::ostringstream oss;
  // oss << std::scientific << max_list_size;

  // std::string scientificString = oss.str();
  // std::ofstream outputFile(outputFilePath + "DLD_" + scientificString + "_" +
  //                          std::to_string(double(snr_dB))+ ".txt");
  // if (!outputFile.is_open()) {
  //   std::cerr << "Failed to open the file for writing." << std::endl;
  // }

  /// Building Codec
  CodeInformation code;
  code.k = 1;
  code.n = 3;
  code.v = V;
  code.list_size = 5;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {AB, AC, DC};  // octal
  ViterbiCodec codec(code);

  CodeInformation code_1;
  // 3325
  // x^10 + x^9 + x^7 + x^6 + x^4 + x^2 + 1
  // CC:(x^5 + x^4 + 1) 110001(binary), 61(octal), 49(decimal)
  // CRC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)

  // 3237
  // x^10 + x^9 + x^7 + x^4 + x^3 + x^2 + x + 1
  // CC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CRC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)
  code_1.k = 1;
  code_1.n = 2;
  code_1.v = HALF_V;
  code_1.list_size = 5e4;  // This is only useful when using list decoding
                           // functions in ViterbiCodec
  code_1.crc_dec = CRC_A;  // 37
  code_1.crc_length = CONSTRAINT_LENGTH;
  code_1.generator_poly = {B, C};  // octal

  CodeInformation code_2;
  // 3237
  // x^10 + x^9 + x^7 + x^4 + x^3 + x^2 + x + 1
  // CRC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)

  // 2711
  // x^10 + x^8 + x^7 + x^6 + x^3 + 1
  // CRC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CC: (x^5 + x^4 + x^2 + x + 1) 110111(binary), 67(octal), 55(decimal)
  code_2.k = 1;
  code_2.n = 2;
  code_2.v = HALF_V;
  code_2.list_size = 5e4;  // This is only useful when using list decoding
                           // functions in ViterbiCodec
  code_2.crc_dec = CRC_C;
  code_2.crc_length = CONSTRAINT_LENGTH;
  code_2.generator_poly = {A, D};

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  std::vector<CodeInformation> dld_codes = {code_1, code_2};
  DualListDecoder DLD(code, dld_codes, max_list_size);

  // Simulation begin
  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 64;

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_list_exceed_vec;
  std::vector<double> DLD_error_vec;

  std::cout << "Running DLD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << snr_dB << std::endl;

  std::cout << "Now working on snr: " << snr_dB << "-------------------"
            << std::endl;

  std::vector<int> expected_list_ranks = {0, 0};

  int DLD_correct = 0;
  int DLD_list_exceeded = 0;
  int DLD_error = 0;

  int SSV_error = 0;

  // vector to keep track of list sizes for both decoders
  std::vector<int> DLD_list_0_size;
  std::vector<int> DLD_list_1_size;

  int number_of_trials = 0;
  int number_of_errors = 0;

  std::vector<double> received_to_decoded_dist;
  std::vector<double> received_to_correct_dist;

  std::vector<double> metric_0((K + code.v) * 3, 0.0);
  std::vector<double> metric_1((K + code.v) * 3, 0.0);

  std::vector<std::chrono::milliseconds> timeDurations(
      4, std::chrono::milliseconds(0));
  std::chrono::milliseconds softDurations(0);

  // while (number_of_errors < max_errors) {
  while (number_of_trials < TRIALS) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
      std::cout << "Trial number: " << number_of_trials << std::endl;
      std::cout << "Current number of errors: " << number_of_errors
                << ", SSV errors: " << SSV_error
                << std::endl;
    }

    std::vector<int> msg;
    for (int i = 0; i < num_bits; ++i) {
      int random_bit = rand() % 2;
      msg.push_back(random_bit);
    }

    std::vector<int> encoded_msg = codec.encodeZTCC(msg);

    assert(encoded_msg.size() == (msg.size() + code.v) * 3);

    std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);

    std::vector<double> received_signal =
        AWGN::addNoise(modulated_signal, snr_dB);

    // if (number_of_trials < 5) {
    //   continue;
    // }
    // if (number_of_trials > 5) {
    //   break;
    // }

    if (PUNCTURE) {
      for (int i = 0; i < received_signal.size(); i += 6) {
        received_signal[i + PUNC_1] = 0;
        received_signal[i + PUNC_3] = 0;
      }
    }

    for (int i = 0; i < received_signal.size(); i++) {
      metric_0[i] = std::pow(1.0 - received_signal[i], 2);
      metric_1[i] = std::pow(-1.0 - received_signal[i], 2);
    }

    // // Using ZTCC list decoding to verify my result
    // std::vector<MessageInformation> output =
    // codec.ZTCCListDecoding_fullInformation_NoConstraint(received_signal);
    // outputFile << "printing single ztcc full information list decoding
    // results: " << std::endl; for (int i = 0; i < output.size(); ++i) {
    //   outputFile << " " << i << "th path  (" << "list rank = " <<
    //   output[i].list_rank << "): "; CodecUtils::outputMat(output[i].message,
    //   outputFile); outputFile << " with pathMetric = " <<
    //   output[i].path_metric; outputFile << std::endl;
    // }
    // outputFile << std::endl;

    // // unleaver to unleave the bits from received_signal
    // std::vector<double> received_codec_2;
    // std::vector<double> received_codec_1;
    // for (size_t i = 0; i < received_signal.size(); ++i) {
    //   if (i % 3 == 0) {
    //     received_codec_2.push_back(received_signal[i]);
    //     received_codec_2.push_back(received_signal[i+1]);
    //   } else if (i % 3 == 1) {
    //     received_codec_1.push_back(received_signal[i]);
    //     received_codec_1.push_back(received_signal[i+1]);
    //   } else {
    //     continue;
    //   }
    // }
    // assert(received_codec_2.size() == 148);
    // assert(received_codec_1.size() == 148);

    // std::vector<MessageInformation> output_1 =
    // codec_1.ZTCCListDecoding_fullInformation_WithConstraint(received_codec_2);
    // outputFile << "printing DLD1 full information list decoding results: " <<
    // std::endl; for (int i = 0; i < output_1.size(); ++i) {
    //   outputFile << " " << output_1[i].crc_passing_rank << "th crc-passing
    //   path  (" << "list rank = " << output_1[i].list_rank << "): ";
    //   CodecUtils::outputMat(output_1[i].message, outputFile);
    //   outputFile << " with pathMetric = " << output_1[i].path_metric;
    //   outputFile << std::endl;
    // }
    // outputFile << std::endl;

    // std::vector<MessageInformation> output_2 =
    // codec_2.ZTCCListDecoding_fullInformation_WithConstraint(received_codec_1);
    // outputFile << "printing DLD2 full information list decoding results: " <<
    // std::endl; for (int i = 0; i < output_2.size(); ++i) {
    //   outputFile << " " << output_2[i].crc_passing_rank << "th crc-passing
    //   path  (" << "list rank = " << output_2[i].list_rank << "): ";
    //   CodecUtils::outputMat(output_2[i].message, outputFile);
    //   outputFile << " with pathMetric = " << output_2[i].path_metric;
    //   outputFile << std::endl;
    // }
    // outputFile << std::endl;

    // DLD DECODING
    DLDInfo output_DLD = DLD.LookAheadDecode_SimpleAlternate_StopOnceMatchFound_WithListSizeExceeded_HalfMetricOnSharedSymbols(
        received_signal, timeDurations);

    // we save the list ranks
    DLD_list_0_size.push_back(output_DLD.list_ranks[0]);
    DLD_list_1_size.push_back(output_DLD.list_ranks[1]);

    if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
      // CASE 1
      // correct decoding
      received_to_correct_dist.push_back(output_DLD.combined_metric);
      received_to_decoded_dist.push_back(output_DLD.combined_metric);
      DLD_correct++;
      // MessageInformation output_SSV =
      //     codec.softViterbiDecoding(received_signal, timeDurations[1]);
      // if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
      //   // std::cout << "trial number: " << number_of_trials << std::endl;
      //   // std::cout << "DLD makes an undectected error, now trying SSV"
      //   //           << std::endl;
      //   std::cout << "SSV decodes correctly" << std::endl;
      //   std::cout << "DLD metric: " << output_DLD.combined_metric << std::endl;
      //   std::cout << "SSV metrics: " << output_SSV.path_metric << std::endl;
      // }
    } else if (output_DLD.message == std::vector<int>(64, -1)) {
      // CASE 2
      // list size exceeded
      DLD_list_exceeded++;
      number_of_errors++;

    } else {
      // CASE 3
      // if there is an error, we save the correct message, the incorrectly
      // decoded message, and the received signal and record the list index.
      DLD_error++;
      number_of_errors++;
      // std::cout << "DLD decoding error!" << std::endl;
      MessageInformation output_SSV =
          codec.softViterbiDecoding(received_signal, timeDurations[1]);
      if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
        // std::cout << "trial number: " << number_of_trials << std::endl;
        // std::cout << "DLD makes an undectected error, now trying SSV"
        //           << std::endl;
        std::cout << "SSV decodes correctly" << std::endl;
        std::cout << "DLD metric: " << output_DLD.combined_metric << std::endl;
        std::cout << "SSV metrics: " << output_SSV.path_metric << std::endl;

        // DualListDecoder DLD_Extra_LS(code, dld_codes, 10*max_list_size);
        // DLDInfo output_DLD_LS = DLD_Extra_LS.LookAheadDecode_SimpleAlternate_StopOnceMatchFound_WithListSizeExceeded_HalfMetricOnSharedSymbols(
        //     received_signal, timeDurations);
        // if (CodecUtils::areVectorsEqual(output_DLD_LS.message, msg)) {
        //   std::cout << "DLD Extra LS metric: "
        //             << output_DLD_LS.combined_metric << ", DLD Extra LS decoding success!" << std::endl;
        // } else {
        //   std::cout << "DLD Extra LS decoding error!" << std::endl;
        // }
      } else {
        SSV_error++;
      }
    }

    // update list ranks
    for (int i = 0; i < expected_list_ranks.size(); ++i) {
      expected_list_ranks[i] += output_DLD.list_ranks[i];
    }
  }


  for (int j = 0; j < expected_list_ranks.size(); ++j) {
    expected_list_ranks[j] /= number_of_trials;
  }
  std::cout << "DLD correct decoding: " << DLD_correct
            << " , Percentage: " << (double)DLD_correct / number_of_trials
            << std::endl;
  std::cout << "DLD Exceeded list:" << DLD_list_exceeded
            << ", Percentage: " << (double)DLD_list_exceeded / number_of_trials
            << std::endl;
  std::cout << "DLD wrong decoding: " << DLD_error
            << " , Percentage: " << (double)DLD_error / number_of_trials
            << std::endl;

  std::cout << "SSV wrong decoding: " << SSV_error
          << " , Percentage: " << (double)SSV_error / number_of_trials
          << std::endl;

  std::cout << "expected list ranks: [ " << expected_list_ranks[0] << ", "
            << expected_list_ranks[1] << " ]" << std::endl;

  // distance spectrum statistics
  // std::cout << "List ranks spectrum data: " << std::endl;
  // std::cout << "  List 0: ";
  // CodecUtils::printMat(DLD_list_0_size);
  // std::cout << std::endl;
  // std::cout << "  List 1: ";
  // CodecUtils::printMat(DLD_list_1_size);
  // std::cout << std::endl;

  DLD_correct_vec.push_back((double)DLD_correct / number_of_trials);
  DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
  DLD_error_vec.push_back((double)DLD_error / number_of_trials);

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
  // outputFile.close();
}

mixed_info mixedDualListExperiment_rate_1_3(double snr_dB,
                                            int maximum_list_size,
                                            int max_errors) {
  /*
  function description here
  TODO:

  */
  // output path
  // std::string outputFilePath = "../output/list_size_exceeded/";
  // std::ostringstream oss;
  // oss << std::scientific << maximum_list_size;

  // std::string scientificString = oss.str();
  // std::ofstream outputFile(outputFilePath + "DLD_" + scientificString + "_" +
  //                          std::to_string(double(snr_dB))+ ".txt");
  // if (!outputFile.is_open()) {
  //   std::cerr << "Failed to open the file for writing." << std::endl;
  // }

  // return
  mixed_info result;

  /// Building Codec
  CodeInformation code;
  code.k = 1;
  code.n = 3;
  code.v = V;
  code.list_size = 5;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {AB, AC, DC};  // octal
  ViterbiCodec codec(code);

  CodeInformation code_1;
  // 3325
  // x^10 + x^9 + x^7 + x^6 + x^4 + x^2 + 1
  // CC:(x^5 + x^4 + 1) 110001(binary), 61(octal), 49(decimal)
  // CRC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)

  // 3237
  // x^10 + x^9 + x^7 + x^4 + x^3 + x^2 + x + 1
  // CC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CRC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)
  code_1.k = 1;
  code_1.n = 2;
  code_1.v = 5;
  code_1.list_size = 5e4;  // This is only useful when using list decoding
                           // functions in ViterbiCodec
  code_1.crc_dec = CRC_A;  // 37
  code_1.crc_length = 6;
  code_1.generator_poly = {B, C};  // octal

  CodeInformation code_2;
  // 3237
  // x^10 + x^9 + x^7 + x^4 + x^3 + x^2 + x + 1
  // CRC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CC: (x^5 + x^2 + 1) 100101(binary), 45(octal), 37(decimal)

  // 2711
  // x^10 + x^8 + x^7 + x^6 + x^3 + 1
  // CRC:(x^5 + x^4 + x + 1) 110011(binary), 63(octal), 51(decimal)
  // CC: (x^5 + x^4 + x^2 + x + 1) 110111(binary), 67(octal), 55(decimal)
  code_2.k = 1;
  code_2.n = 2;
  code_2.v = 5;
  code_2.list_size = 5e4;  // This is only useful when using list decoding
                           // functions in ViterbiCodec
  code_2.crc_dec = CRC_C;
  code_2.crc_length = 6;
  code_2.generator_poly = {A, D};

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  std::vector<CodeInformation> dld_codes = {code_1, code_2};
  DualListDecoder DLD(code, dld_codes, maximum_list_size);

  // Simulation begin
  srand(seed);
  noise_gen.seed(82);

  std::vector<double> DLD_correct_vec;
  std::vector<double> DLD_list_exceed_vec;
  std::vector<double> DLD_error_vec;

  std::cout << "Running DLD tests for v = " << code.v
            << " test, with generator poly: " << code.generator_poly[0] << ", "
            << code.generator_poly[1] << " SNR = " << snr_dB << std::endl;

  std::cout << "DLD order: g1 = " << code_1.generator_poly[0]
            << ", g2 = " << code_2.generator_poly[0] << std::endl;

  std::vector<std::chrono::milliseconds> stepTimeDurations(
      4, std::chrono::milliseconds(0));
  std::chrono::milliseconds overall_decoding_time_debug(0);

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

  std::vector<double> metric_0((K + code.v) * 3, 0.0);
  std::vector<double> metric_1((K + code.v) * 3, 0.0);

  // measure the time of each operation
  // step1: SSV_list_1: time taken to add-compare-select and perform the initial
  // trace-back step2: C_trace_1: time taken for additional traaceback
  // opearations required by SLVD step3: C_insert_1: time taken of inserting new
  // elements to maintain an ordered list of path metric differences

  // std::chrono::milliseconds softDurations(0);

  // while (number_of_errors < max_errors) {
  Stopwatch overall_decoding_time_debug_sw;
  overall_decoding_time_debug_sw.tic();
  
  while (number_of_trials < TRIALS) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
      std::cout << "Trial number: " << number_of_trials << std::endl;
      std::cout << "Current number of errors: " << number_of_errors
                << std::endl;
    }

    std::vector<int> msg;
    for (int i = 0; i < K; ++i) {
      int random_bit = rand() % 2;
      msg.push_back(random_bit);
    }

    // coding
    std::vector<int> encoded_msg = codec.encodeZTCC(msg);

    assert(encoded_msg.size() == (msg.size() + code.v) * 3);

    std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);

    std::vector<double> received_signal =
        AWGN::addNoise(modulated_signal, snr_dB);

    if (PUNCTURE) {
      for (int i = 0; i < received_signal.size(); i += 6) {
        received_signal[i + PUNC_1] = 0;
        received_signal[i + PUNC_3] = 0;
      }
    }

    for (int i = 0; i < received_signal.size(); i++) {
      metric_0[i] = std::pow(1.0 - received_signal[i], 2);
      metric_1[i] = std::pow(-1.0 - received_signal[i], 2);
    }

    // continue;
    

    // DLD DECODING
    DLDInfo output_DLD = DLD.LookAheadDecode_SimpleAlternate_rate_1_2(
        received_signal, stepTimeDurations, metric_0, metric_1);



    // if DLD declares List size exceeded
    // then resort to soft viterbi decoding
    if (output_DLD.message == std::vector<int>(64, -1)) {

      // std::cout << "List size EXCEEDED! Running SSV Nows" << " " << SSV_correct << std::endl;
      MessageInformation output_SSV =
          codec.softViterbiDecoding(received_signal, stepTimeDurations[1]);
      if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
        SSV_correct++;
      } else {
        SSV_error++;
        number_of_errors++;
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

      // MessageInformation output_SSV =
      //     codec.softViterbiDecoding(received_signal, stepTimeDurations[1]);
      // if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
      //   std::cout << "DLD makes an undectected error, now trying SSV"
      //             << std::endl;
      //   std::cout << "SSV decodes correctly" << std::endl;
      //   std::cout << "SSV metrics: " << output_SSV.path_metric << std::endl;
      //   std::cout << number_of_trials << std::endl;
      // }
    }

    // update list ranks
    for (int i = 0; i < expected_list_ranks.size(); ++i) {
      expected_list_ranks[i] += output_DLD.list_ranks[i];
    }


  }

  overall_decoding_time_debug_sw.toc();
  overall_decoding_time_debug += overall_decoding_time_debug_sw.getElapsed();
  overall_decoding_time_debug_sw.reset();

  ////////  DLD Decoding Statistics /////////
  // debug
  std::cout << "overall decoding time debug: " << overall_decoding_time_debug.count() << " milliseconds" << std::endl;

  // record the time taken by each step
  for (size_t i = 0; i < stepTimeDurations.size(); ++i) {
    std::cout << "Step " << i + 1
              << " Runtime: " << stepTimeDurations[i].count() << " milliseconds"
              << std::endl;
  }

  std::cout << "Out of " << number_of_trials << " experiments, "
            << DLD_list_0_size.size()
            << " experiments were recorded to have DLD list sizes."
            << std::endl;

  for (int j = 0; j < expected_list_ranks.size(); ++j) {
    expected_list_ranks[j] /= DLD_list_0_size.size();
  }
  std::cout << "DLD correct decoding: " << DLD_correct
            << " , Percentage: " << (double)DLD_correct / number_of_trials
            << std::endl;
  std::cout << "DLD Exceeded list:" << DLD_list_exceeded
            << ", Percentage: " << (double)DLD_list_exceeded / number_of_trials
            << std::endl;
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

  DLD_correct_vec.push_back((double)(DLD_correct + SSV_correct) /
                            number_of_trials);
  DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
  DLD_error_vec.push_back((double)DLD_error + SSV_error / number_of_trials);

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
  result.DLD_decoded_portion =
      (double)DLD_list_0_size.size() / number_of_trials;
  result.DLD_expected_list_sizes = expected_list_ranks;
  return result;
}

std::chrono::milliseconds softViterbiExperiment_rate_1_3(double snr_dB,
                                                         int max_errors) {
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
  code.n = 3;
  code.v = 10;
  code.list_size = 5;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {AB, AC, DC};  // octal
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

  // while (number_of_errors < max_errors) {
  while (number_of_trials < TRIALS) {
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

    assert(encoded_msg.size() == (msg.size() + code.v) * 3);

    std::vector<int> modulated_signal = BPSK::modulate(encoded_msg);
    // outputFile << "Printing modulated signal: " << std::endl;
    // CodecUtils::outputMat(modulated_signal, outputFile);
    // outputFile << std::endl;

    std::vector<double> received_signal =
        AWGN::addNoise(modulated_signal, snr_dB);
    // outputFile << "Printing received signal: " << std::endl;
    // CodecUtils::outputMat(received_signal, outputFile);
    // outputFile << std::endl;

    if (PUNCTURE) {
      for (int i = 0; i < received_signal.size(); i += 6) {
        received_signal[i + PUNC_1] = 0;
        received_signal[i + PUNC_3] = 0;
      }
    }

    MessageInformation output_ML =
        codec.softViterbiDecoding(received_signal, softDurations);

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
            << " , Percentage: " << (double)Correct_decoding / number_of_trials
            << std::endl;
  std::cout << "std viterbi wrong decoding: " << ML_decoding_error
            << " , Percentage: " << (double)ML_decoding_error / number_of_trials
            << std::endl;
  correct_decoding_snr.push_back((double)Correct_decoding / number_of_trials);
  ML_decoding_error_snr.push_back((double)ML_decoding_error / number_of_trials);

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