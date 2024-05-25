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
#include <tuple>

#include "../include/CONSTANTS.h"
#include "../include/dualListDecoder.h"
#include "../include/trialListDecoder.h"
#include "../include/viterbiDecoder.h"
#include "../include/stopWatch.h"
#include "../include/viterbiCodec.h"
// #include "mat.h"

static std::random_device rd{};
int seed = 33;                      // 47
static std::mt19937 msg_gen(seed);



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

std::vector<double> ComputeSquaredDifferences(
    const std::vector<int>& vector1, const std::vector<double>& vector2) {
  // Check if the vectors have the same size
  if (vector1.size() != vector2.size()) {
    // You can handle this error in your preferred way, e.g., throw an exception
    throw std::invalid_argument("Vectors must have the same size");
  }

  // Calculate squared differences
  std::vector<double> squaredDifferences;
  squaredDifferences.reserve(vector1.size());  // Reserve space for efficiency

  for (std::size_t i = 0; i < vector1.size(); ++i) {
    double diff = vector1[i] - vector2[i];
    squaredDifferences.push_back(diff * diff);
  }

  return squaredDifferences;
}

double SumGroupIndexElements(const std::vector<double>& inputVector,
                             std::size_t groupLength, std::size_t groupIndex) {
  // Check if the length of the vector is a multiple of group length

  if (inputVector.size() % groupLength != 0) {
    throw std::invalid_argument(
        "Vector size must be a multiple of group length");
  }

  // Sum together every element of the specified group index in every group
  double sum = 0.0;
  for (std::size_t i = groupIndex; i < inputVector.size(); i += groupLength) {
    sum += inputVector[i];
  }

  return sum;
}
}  // namespace


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

/*
Functions that end with Experiment is a stand-alone experiment designed to run
by itself Functions taht end with Decoding is a decoding function that takes in
a single received message and
*/
std::chrono::milliseconds softViterbiExperiment(double snr_dB, int max_errors);
void dualListExperiment_rate_1_3(double snr_dB, int max_list_size,
                                 int max_errors);
mixed_info mixedDualListExperiment_rate_1_3(double snr_dB,
                                            int maximum_list_size,
                                            int max_errors);
std::chrono::milliseconds softViterbiExperiment_rate_1_3(double snr_dB,
                                                         int max_errors);
void TimeAndComplexitySimulation_rate_1_3(double SNR);
void TimeAndComplexitySimulation_SSV_DLD(double SNR);
std::vector<std::chrono::milliseconds> dualListDecodeTime(double snr_dB, int max_list_size,
                                 int max_errors);

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

  ////////////////////// DSU/Soft Viterbi Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   softViterbiExperiment(snr_dB, 5);
  // }

  ////////////////////// RATE 1/3 Dual List Decoding Experiment
  /////////////////////////
  for (double snr_dB : SNR_dB) {
    dualListExperiment_rate_1_3(snr_dB, MAX_LIST_SIZE, 20);
  }

  ////////////////////// RATE 1/3 DSU/Soft Viterbi Experiment
  /////////////////////////
  // for (double snr_dB : SNR_dB) {
  //   softViterbiExperiment_rate_1_3(snr_dB, 5);
  // }

  ////////////////////// Time and Complexity Experiment //////////////////////
  // for (double snr_dB : SNR_dB) {
  //   TimeAndComplexitySimulation_rate_1_3(snr_dB);
  // }

  // for (double snr_dB : SNR_dB) {
  //   TimeAndComplexitySimulation_SSV_DLD(snr_dB);
  // }

  ////////////////////// Find D_free //////////////////////
  // FindDFree_ADFree();
  return 0;
}

void TimeAndComplexitySimulation_SSV_DLD(double SNR) {
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
  std::vector<std::chrono::milliseconds> DLD_overall_time;

  std::vector<std::chrono::milliseconds> standard_soft_viterbi_time;
  for (int DLD_maximum_list_size : max_list_size_vector) {
    std::cout << "Running experiments for max list size  = "
              << DLD_maximum_list_size << std::endl;
    std::vector<std::chrono::milliseconds> result = dualListDecodeTime(
        SNR, DLD_maximum_list_size, max_errors);
    // record the results
    mixed_ssv_time.push_back(result[0]);
    reconstruct_trellis_time.push_back(result[1]);
    mixed_insertion_time.push_back(result[2]);
    DLD_overall_time.push_back(result[3]);

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

  std::cout << "printing DLD overall time: [";
  for (std::chrono::milliseconds time : DLD_overall_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;

  std::cout << "printing standard soft viterbi time: [";
  for (std::chrono::milliseconds time : standard_soft_viterbi_time) {
    std::cout << time.count() << ", ";
  }
  std::cout << "]" << std::endl;
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

void dualListExperiment_rate_1_3(double snr_dB, int max_list_size,
                                 int max_errors) {

  std::string outputFilePath = "../output/single_symbol_metric/";
  std::ostringstream oss;
  oss << std::scientific << max_list_size;
  std::string scientificString = oss.str();
  std::ofstream outputFile(outputFilePath + "DLD_" + scientificString + "_" +
                           std::to_string(double(snr_dB)) + ".txt");
  if (!outputFile.is_open()) {
    std::cerr << "Failed to open the file for writing." << std::endl;
  }

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
  code_1.k = 1;
  code_1.n = 2;
  code_1.v = HALF_V;
  code_1.list_size = 5e4;
  code_1.crc_dec = CRC_A;
  code_1.crc_length = CONSTRAINT_LENGTH;
  code_1.generator_poly = {B, C};

  CodeInformation code_2;
  code_2.k = 1;
  code_2.n = 2;
  code_2.v = HALF_V;
  code_2.list_size = 5e4;
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
  int DLD_half_on_shared_symbol_error = 0;
  int DLD_double_match_error = 0;
  int SSV_error = 0;
  int number_of_trials = 0;
  int number_of_errors = 0;

  std::vector<std::tuple<double, double, double>> Given_SSV_succeeded_DLD_failed_symbol_metrics;
  std::vector<std::tuple<double, double, double>> Given_SSV_succeeded_SSV_success_symbol_metrics;
  std::vector<std::tuple<double, double, double>> DLD_succeeded_symbol_metrics; 

  // std::vector<int> DLD_full_on_shared_metric_expected_list_ranks = {0, 0};
  // std::vector<int> DLD_half_on_shared_metric_expected_list_ranks = {0, 0};

  std::vector<std::chrono::milliseconds> timeDurations(
      4, std::chrono::milliseconds(0));
  std::chrono::milliseconds softDurations(0);

  // while (number_of_errors < max_errors) {
  while (number_of_trials < TRIALS) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
      std::cout << "Trial number: " << number_of_trials << std::endl;
      std::cout << "Current number of errors: " << number_of_errors
                << ", SSV errors: " << SSV_error << std::endl;
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

    if (PUNCTURE) {
      for (int i = 0; i < received_signal.size(); i += 6) {
        received_signal[i + PUNC_1] = 0;
        received_signal[i + PUNC_3] = 0;
      }
    }

    // DLD DECODING
    DLDInfo output_DLD =
        DLD.DLD_BAM(
            received_signal, timeDurations);

    if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
      // CASE 1
      // correct decoding

      DLD_correct++;

      DLD_succeeded_symbol_metrics.push_back(output_DLD.symbol_metrics);
      
      /////////////////////////////////////////////////  DLD on shared metric /////////////////////////////////////////////////
      
      // std::cout << "shared metric not addressed List size symbol metric: "
      //           << std::get<0>(output_DLD.symbol_metrics) << ", "
      //           << std::get<1>(output_DLD.symbol_metrics) << ", "
      //           << std::get<2>(output_DLD.symbol_metrics) << std::endl;

      // DLD_full_on_shared_metric_expected_list_ranks[0] += output_DLD.list_ranks[0];
      // DLD_full_on_shared_metric_expected_list_ranks[1] += output_DLD.list_ranks[1];

      // DLDInfo output_DLD_half_shared_metric = DLD.DLD_BAM_Half_Metric_on_Shared(
      //                                             received_signal, timeDurations);
      // if (CodecUtils::areVectorsEqual(output_DLD_half_shared_metric.message, msg)) {
      //   // std::cout << "shared metric addressed List size symbol metric: "
      //   //           << std::get<0>(output_DLD_half_shared_metric.symbol_metrics) << ", "
      //   //           << std::get<1>(output_DLD_half_shared_metric.symbol_metrics) << ", "
      //   //           << std::get<2>(output_DLD_half_shared_metric.symbol_metrics) << std::endl;
        
      //   DLD_half_on_shared_metric_expected_list_ranks[0] += output_DLD_half_shared_metric.list_ranks[0];
      //   DLD_half_on_shared_metric_expected_list_ranks[1] += output_DLD_half_shared_metric.list_ranks[1];
        
      // } else {
      //   DLD_half_on_shared_symbol_error++;
      //   std::cout << "shared symbol method decode incorrectly." << std::endl;
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
      // std::cout << "1x List size symbol metric: "
      //           << std::get<0>(output_DLD.symbol_metrics) << ", "
      //           << std::get<1>(output_DLD.symbol_metrics) << ", "
      //           << std::get<2>(output_DLD.symbol_metrics) << std::endl;
      // std::cout << "DLD list 1 metric: "
      //           << std::get<0>(output_DLD.symbol_metrics) +
      //                  std::get<1>(output_DLD.symbol_metrics)
      //           << ", DLD List 2 metric:"
      //           << std::get<1>(output_DLD.symbol_metrics) +
      //                  std::get<2>(output_DLD.symbol_metrics)
      //           << std::endl;

      // DLD_full_on_shared_metric_expected_list_ranks[0] += output_DLD.list_ranks[0];
      // DLD_full_on_shared_metric_expected_list_ranks[1] += output_DLD.list_ranks[1];

      /////////////////////////////////////////////////  DLD on shared metric /////////////////////////////////////////////////

      // DLDInfo output_DLD_half_shared_metric = DLD.DLD_BAM_Half_Metric_on_Shared(
      //                                               received_signal, timeDurations);
      //   if (CodecUtils::areVectorsEqual(output_DLD_half_shared_metric.message, msg)) {
      //     // std::cout << "shared metric addressed List size symbol metric: "
      //     //           << std::get<0>(output_DLD_half_shared_metric.symbol_metrics) << ", "
      //     //           << std::get<1>(output_DLD_half_shared_metric.symbol_metrics) << ", "
      //     //           << std::get<2>(output_DLD_half_shared_metric.symbol_metrics) << std::endl;
          
      //     DLD_half_on_shared_metric_expected_list_ranks[0] += output_DLD_half_shared_metric.list_ranks[0];
      //     DLD_half_on_shared_metric_expected_list_ranks[1] += output_DLD_half_shared_metric.list_ranks[1];
      //     std::cout << "shared symbol method decode correctly." << std::endl;
      //   } else {
      //     DLD_half_on_shared_symbol_error++;
      //     std::cout << "shared symbol method decode incorrectly." << std::endl;
      //   }

      //////////////////////////////////////////  SSV /////////////////////////////////////////////////

      // MessageInformation output_SSV =
      //     codec.softViterbiDecoding(received_signal, timeDurations[1]);
      // if (CodecUtils::areVectorsEqual(output_SSV.message, msg)) {
      //   // std::cout << "SSV decodes correctly" << std::endl;

      //   std::vector<int> reencoded_msg = codec.encodeZTCC(output_SSV.message);
      //   assert(reencoded_msg.size() == received_signal.size());

      //   std::vector<int> reencoded_signal = BPSK::modulate(reencoded_msg);

      //   std::vector<double> reencode_squared_differences =
      //       ComputeSquaredDifferences(reencoded_signal, received_signal);

      //   double symbol1_metric =
      //       SumGroupIndexElements(reencode_squared_differences, 3, 0);
      //   double symbol2_metric =
      //       SumGroupIndexElements(reencode_squared_differences, 3, 1);
      //   double symbol3_metric =
      //       SumGroupIndexElements(reencode_squared_differences, 3, 2);
      //   auto symbol_metrics =
      //       std::make_tuple(symbol1_metric, symbol2_metric, symbol3_metric);
      //   Given_SSV_succeeded_SSV_success_symbol_metrics.push_back(symbol_metrics);
      //   Given_SSV_succeeded_DLD_failed_symbol_metrics.push_back(output_DLD.symbol_metrics);
      // } else {
      //   // std::cout << "SSV decodes incorrectly" << std::endl;
      //   SSV_error++;
      // }

      ////////////////////////////////////////// DLD Double Match /////////////////////////////////////////////////

      DLDInfo output_DLD_double_match =
        DLD.DLD_BAM_double_match(
            received_signal, timeDurations);
      if (output_DLD_double_match.message == msg) {
        std::cout << "DLD double match correct decoding." << std::endl;
      } else {
        std::cout << "DLD double match incorrect decoding." << std::endl;
        DLD_double_match_error++;
      }

    }
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
  
  // std::cout << "DLD half on shared symbol error: " << DLD_half_on_shared_symbol_error
  //           << " , Percentage: " << (double)DLD_half_on_shared_symbol_error / number_of_trials
  //           << std::endl;

  std::cout << "DLD double match error: " << DLD_double_match_error
            << " , Percentage: " << (double)DLD_double_match_error / number_of_trials
            << std::endl;

  std::cout << "SSV wrong decoding: " << SSV_error
            << " , Percentage: " << (double)SSV_error / number_of_trials
            << std::endl;

  // DLD_correct_vec.push_back((double)DLD_correct / number_of_trials);
  // DLD_list_exceed_vec.push_back((double)DLD_list_exceeded / number_of_trials);
  // DLD_error_vec.push_back((double)DLD_error / number_of_trials);

  //  /////////////////////////////////////////////////  DLD single symbol metric /////////////////////////////////////////////////
  // if (outputFile.is_open()) {
  //   outputFile << "DLD Failed, the incorrect DLD metrics:" << std::endl;
  //   for (const auto& metrics : Given_SSV_succeeded_DLD_failed_symbol_metrics) {
  //     outputFile << std::get<1>(metrics) << std::endl;
  //   }
  //   outputFile << std::endl;

  //   outputFile << "DLD Failed, the ML SSV metrics:" << std::endl;
  //   for (const auto& metrics : Given_SSV_succeeded_SSV_success_symbol_metrics) {
  //     outputFile << std::get<1>(metrics) << std::endl;
  //   }
  //   outputFile << std::endl;

  //   outputFile << "DLD Succeeded Symbol Metrics:" << std::endl;
  //   for (const auto& metrics : DLD_succeeded_symbol_metrics) {
  //     outputFile << std::get<1>(metrics) << std::endl;
  //   }
  //   outputFile << std::endl;
  // } else {
  //   std::cout << "Failed to open output file." << std::endl;
  // }


   /////////////////////////////////////////////////  DLD on shared metric /////////////////////////////////////////////////
  // for (int& value : DLD_full_on_shared_metric_expected_list_ranks) {
  //   value = value / DLD_correct;
  // }
  // std::cout << "DLD full on shared metric expected list ranks: "
  //           << DLD_full_on_shared_metric_expected_list_ranks[0] << ", "
  //           << DLD_full_on_shared_metric_expected_list_ranks[1] << std::endl;
  
  // for (int& value : DLD_half_on_shared_metric_expected_list_ranks) {
  //   value = value / DLD_correct;
  // }
  // std::cout << "DLD half on shared metric expected list ranks: "
  //           << DLD_half_on_shared_metric_expected_list_ranks[0] << ", "
  //           << DLD_half_on_shared_metric_expected_list_ranks[1] << std::endl;
  

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


std::vector<std::chrono::milliseconds> dualListDecodeTime(double snr_dB, int max_list_size,
                                 int max_errors) {

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
  code_1.k = 1;
  code_1.n = 2;
  code_1.v = HALF_V;
  code_1.list_size = 5e4;
  code_1.crc_dec = CRC_A;
  code_1.crc_length = CONSTRAINT_LENGTH;
  code_1.generator_poly = {B, C};

  CodeInformation code_2;
  code_2.k = 1;
  code_2.n = 2;
  code_2.v = HALF_V;
  code_2.list_size = 5e4;
  code_2.crc_dec = CRC_C;
  code_2.crc_length = CONSTRAINT_LENGTH;
  code_2.generator_poly = {A, D};
  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);
  std::vector<CodeInformation> dld_codes = {code_1, code_2};
  DualListDecoder DLD(code, dld_codes, max_list_size);

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
  std::cout << "Now working on snr: " << snr_dB << "-------------------"
            << std::endl;

  std::vector<int> expected_list_ranks = {0, 0};

  int DLD_correct = 0;
  int DLD_list_exceeded = 0;
  int DLD_error = 0;
  int SSV_error = 0;
  int number_of_trials = 0;
  int number_of_errors = 0;

  std::vector<std::chrono::milliseconds> timeDurations(
      4, std::chrono::milliseconds(0));

  // while (number_of_errors < max_errors) {
  while (number_of_trials < TRIALS) {
    number_of_trials++;

    if (number_of_trials % 2000 == 0) {
      std::cout << "Trial number: " << number_of_trials << std::endl;
      std::cout << "Current number of errors: " << number_of_errors
                << ", SSV errors: " << SSV_error << std::endl;
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

    if (PUNCTURE) {
      for (int i = 0; i < received_signal.size(); i += 6) {
        received_signal[i + PUNC_1] = 0;
        received_signal[i + PUNC_3] = 0;
      }
    }

    // DLD DECODING
    Stopwatch DLD_swatch;
    DLD_swatch.tic();
    DLDInfo output_DLD =
        DLD.DLD_BAM(
            received_signal, timeDurations);
    DLD_swatch.toc();
    timeDurations[3] += DLD_swatch.getElapsed();
    DLD_swatch.reset();


    if (CodecUtils::areVectorsEqual(output_DLD.message, msg)) {
      // CASE 1
      // correct decoding
      DLD_correct++;
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
    }
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

  return timeDurations;
}


// void TrialListDecoderSimulation(double SNR, int max_list_size) {
//   CodeInformation code;
//   code.k = 1;
//   code.n = 3;
//   code.v = V;
//   code.list_size = 5;
//   code.crc_dec = -1;
//   code.crc_length = -1;
//   code.generator_poly = {AB, AC, BC};  // octal
//   ViterbiCodec codec(code);

//   CodeInformation code_1; // AB, AC
//   code_1.k = 1;
//   code_1.n = 2;
//   code_1.v = HALF_V;
//   code_1.list_size = 5e4;
//   code_1.crc_dec = CRC_A;
//   code_1.crc_length = CONSTRAINT_LENGTH;
//   code_1.generator_poly = {B, C};

//   CodeInformation code_2; // AC, BC
//   code_2.k = 1;
//   code_2.n = 2;
//   code_2.v = HALF_V;
//   code_2.list_size = 5e4;
//   code_2.crc_dec = CRC_C;
//   code_2.crc_length = CONSTRAINT_LENGTH;
//   code_2.generator_poly = {A, B};

//   CodeInformation code_3; // AB, BC
//   code_3.k = 1;
//   code_3.n = 2;
//   code_3.v = HALF_V;
//   code_3.list_size = 5e4;
//   code_3.crc_dec = CRC_B;
//   code_3.crc_length = CONSTRAINT_LENGTH;
//   code_3.generator_poly = {A, C};
//   ViterbiCodec codec_1(code_1);
//   ViterbiCodec codec_2(code_2);
//   ViterbiCodec codec_3(code_3);
//   std::vector<CodeInformation> tld_codes = {code_1, code_2, code_3};
//   TrialListDecoder TLD(code, tld_codes, max_list_size);


  
// }