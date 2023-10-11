#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <random>

#include "../include/viterbiCodec.h"
#include "mat.h"

namespace AWGN {

std::vector<double> addNoise(std::vector<int> modulated_signal, double SNR) {
  std::random_device rd;
  std::mt19937 noise_gen(47);  // 82
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
void writeToMat(const std::vector<std::vector<T>>& data, const char* filePath,
                const char* varName) {
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

}  // namespace MatlabUtils

namespace DualListDecoder {

struct DLDInfo {
  double combined_metric;
  std::vector<int> message;
  std::vector<int> list_ranks;
};

// Define a custom comparison function for the priority queue
struct CompareCombinedMetric {
  bool operator()(const DLDInfo& a, const DLDInfo& b) const {
    return a.combined_metric >
           b.combined_metric;  // Lower combined_metric at the top
  }
};

// Function to combine two vectors of MessageInformation and create a priority
// queue of DLDInfo sorted in ascending order of combined_metric
std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric>
combine_maps(const std::vector<MessageInformation>& vec1,
             const std::vector<MessageInformation>& vec2) {
  std::map<std::vector<int>, MessageInformation> map1;
  std::map<std::vector<int>, MessageInformation> map2;

  // Populate map1 with elements from vec1
  for (const MessageInformation& msg : vec1) {
    map1[msg.message] = msg;
  }

  // Populate map2 with elements from vec2
  for (const MessageInformation& msg : vec2) {
    map2[msg.message] = msg;
  }

  std::priority_queue<DLDInfo, std::vector<DLDInfo>, CompareCombinedMetric>
      result_queue;

  // Iterate through the keys in map1
  for (const auto& kvp : map1) {
    // Try to find the same key in map2
    auto it = map2.find(kvp.first);
    if (it != map2.end()) {
      // If found, calculate the sum of path_metrics and record list ranks
      double combined_metric = kvp.second.path_metric + it->second.path_metric;
      DLDInfo dld_info;
      dld_info.combined_metric = combined_metric;
      dld_info.message = kvp.first;
      dld_info.list_ranks = {kvp.second.list_rank, it->second.list_rank};
      result_queue.push(dld_info);
    }
  }

  return result_queue;
}

}  // namespace DualListDecoder

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
  int list_size = 500;

  CodeInformation code;
  code.k = 1;
  code.n = 2;
  code.v = 14;
  code.list_size = 1;
  code.crc_dec = -1;
  code.crc_length = -1;
  code.generator_poly = {56721, 61713};
  ViterbiCodec codec(code);

  CodeInformation code_1;
  // x^14+x^13+x^9+x^8+x^7+x^6+x^3+x+1 =
  //                        CRC: (x^3+x^2+1)
  //                        generator: (x^11+x^8+x^7+x^3+x^2+x+1)
  code_1.k = 1;
  code_1.n = 1;
  code_1.v = 11;
  code_1.list_size = list_size;
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
  code_2.list_size = list_size;
  code_2.crc_dec = 7;
  code_2.crc_length = 3;
  code_2.generator_poly = {17453};

  ViterbiCodec codec_1(code_1);
  ViterbiCodec codec_2(code_2);

  int seed = 47;
  std::mt19937 msg_gen(seed);
  int num_bits = 64;
  
  std::vector<int> expected_list_rank_1;
  std::vector<int> expected_list_rank_2;

  for (double snr_dB : SNR_dB) {
    outputFile << "Now working on snr: " << snr_dB << "-------------------"
          << std::endl;
    std::vector<int> cumulative_list_ranks = {0, 0};
    std::vector<double> avg_list_ranks = {0.0, 0.0};
    int NACK_Errors = 0;
    int Correct_decoding = 0;
    int Undetected_Errors = 0;

    for (int trial = 0; trial < mc_N; ++trial) {

      outputFile << "Trial number: " << trial << std::endl;

      std::vector<int> msg;
      for (int i = 0; i < num_bits; ++i) {
        int random_bit = msg_gen() % 2;
        msg.push_back(random_bit);
      }

      // outputFile << "Printing original message: " << std::endl;
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

      std::vector<double> received_codec_2;
      std::vector<double> received_codec_1;

      // unleaver to unleave the bits from received_signal
      for (size_t i = 0; i < received_signal.size(); ++i) {
        if (i % 2 == 0) {
          received_codec_2.push_back(received_signal[i]);
        } else {
          received_codec_1.push_back(received_signal[i]);
        }
      }

      // outputFile << "For codec 1: ";
      // CodecUtils::outputMat(received_codec_1, outputFile);
      // outputFile << std::endl;

      // outputFile << "For codec 2: ";
      // CodecUtils::outputMat(received_codec_2, outputFile);
      // outputFile << std::endl;

      // Decoding starts
      // outputFile << "Decoding begins -----------------" << std::endl;

      std::vector<MessageInformation> output_1 =
          codec_1.ZTCCListViterbiDecoding(received_codec_1);
      std::vector<MessageInformation> output_2 =
          codec_2.ZTCCListViterbiDecoding(received_codec_2);

      std::vector<MessageInformation> output_3 =
          codec.unconstraintZTCCDecoding(received_signal);

      MessageInformation output_soft = codec.softViterbiDecode(received_signal);

      // outputFile << std::endl;

      // outputFile << "Printing Duo list decoder results: " << std::endl;
      // for (int i = 0; i < output_1.size(); ++i) {
      //   outputFile << " " << i << "th message  ("
      //              << "list rank = " << (output_1[i].list_rank)+1 << "): ";
      //   CodecUtils::outputMat(output_1[i].message, outputFile);
      //   outputFile << " with pathMetric = " << output_1[i].path_metric;
      //   outputFile << "     "
      //              << "(list rank = " << (output_2[i].list_rank)+1 << "): ";
      //   CodecUtils::outputMat(output_2[i].message, outputFile);
      //   outputFile << " with pathMetric = " << output_2[i].path_metric;
      //   outputFile << std::endl;
      // }

      // outputFile << std::endl;

      // outputFile << "printing single 2^14 states list decoding results: "
      //            << std::endl;
      // for (int i = 0; i < output_3.size(); ++i) {
      //   outputFile << " " << i << "th message  ("
      //              << "list rank = " << output_3[i].list_rank << "): ";
      //   CodecUtils::outputMat(output_3[i].message, outputFile);
      //   outputFile << " with pathMetric = " << output_3[i].path_metric;
      //   outputFile << std::endl;
      // }

      // outputFile << "printing 2^14 soft viterbi decoding results: " << std::endl;
      // CodecUtils::outputMat(output_soft.message, outputFile);

      // outputFile << std::endl;
      // outputFile << std::endl;

      // Dual List Decoder
      std::priority_queue<DualListDecoder::DLDInfo,
                          std::vector<DualListDecoder::DLDInfo>,
                          DualListDecoder::CompareCombinedMetric>
          result_queue = DualListDecoder::combine_maps(output_1, output_2);

      // Print the result
      if (result_queue.empty()) {
        // outputFile << "No agreed message" << std::endl;
        NACK_Errors++;
      }
      // outputFile << "number of agreed messages: " << result_queue.size() << std::endl;

      while (!result_queue.empty()) {
        DualListDecoder::DLDInfo dld = result_queue.top();
        result_queue.pop();
        double combined_path_metric = dld.combined_metric;
        std::vector<int> list_ranks = dld.list_ranks;
        std::vector<int> message = dld.message;
        // outputFile << "Agreed message: ";
        // CodecUtils::outputMat(message, outputFile);
        // outputFile << "   ";
        // outputFile << "Combined Path Metric: " << combined_path_metric << ",  ";
        // outputFile << "List Ranks: ";
        for (int r = 0; r < list_ranks.size(); ++r) {
          outputFile << list_ranks[r]+1 << ", ";
          cumulative_list_ranks[r] += list_ranks[r]+1;
        }
        // // check decoding result with soft viterbi decoding
        // outputFile << std::endl;
        // outputFile << "Printing first one to compare: ";
        // CodecUtils::outputMat(message, outputFile);
        // outputFile << std::endl;
        // outputFile << "Printing second one to compare: ";
        // CodecUtils::outputMat(output_soft.message, outputFile);
        // outputFile << std::endl;
        // outputFile << std::endl;

        if (CodecUtils::areVectorsEqual(message, output_soft.message)) {
          // outputFile << "Correct decoding!" << std::endl;
          Correct_decoding++;
        } else {
          // outputFile << "Wrong decoding! Incompatible with soft viterbi decoding" << std::endl;
          Undetected_Errors++;
        }
        
        // outputFile << std::endl;
      }

      // outputFile << std::endl;




      // outputFile << std::endl << std::endl;

      // outputFile << "measuring euclidean distance: " <<
      // CodecUtils::euclideanDistance(received_signal, modulated_signal)<<
      // std::endl;

      // std::vector<int> demodulated_signal =
      // BPSK::demodulate(received_signal); outputFile << "Printing demodulated
      // signal: " << std::endl; CodecUtils::outputMat(demodulated_signal,
      // outputFile)  << " with size: " << demodulated_signal.size() <<
      // std::endl; outputFile << std::endl;

      // std::vector<int> decoded_msg =
      // codec.softViterbiDecode(received_signal).message; outputFile <<
      // "Printing soft decoded message: " << std::endl;
      // CodecUtils::outputMat(decoded_msg, outputFile)  << " with size: " <<
      // decoded_msg.size() << std::endl; outputFile << std::endl;

      // std::vector<int> hard_decoded_msg =
      // codec.viterbiDecode(demodulated_signal).message; outputFile <<
      // "Printing hard decoded message: " << std::endl;
      // CodecUtils::outputMat(hard_decoded_msg, outputFile);
      // outputFile << std::endl;

      // std::vector<MessageInformation> output =
      // codec.listViterbiDecoding(received_signal); outputFile << "Printing
      // list decoder results: " << std::endl; for (int i = 0; i <
      // output.size(); ++i) {
      //   outputFile << output[i].message.size() << " " << i << "th message: ";
      //   CodecUtils::output(output[i].message, outputFile);
      //   outputFile << std::endl;
      // }

      // // truncate the flushing bits
      // decoded_msg.resize(msg.size());
      // hard_decoded_msg.resize(msg.size());
      // outputFile << "Measuring hamming distance: " <<
      // CodecUtils::hammingDistance(msg, hard_decoded_msg) << std::endl;
    }
    // outputFile << std::endl << std::endl;
    // outputFile << "Average list ranks for list 1 and 2: ";
    for (int i = 0; i < avg_list_ranks.size(); ++i) {
      avg_list_ranks[i] = (double) cumulative_list_ranks[i] / mc_N;
    }
    // CodecUtils::outputMat(avg_list_ranks, outputFile);

    // outputFile << " NACK Errors: " << NACK_Errors << std::endl;

    expected_list_rank_1.push_back(avg_list_ranks[0]);
    expected_list_rank_2.push_back(avg_list_ranks[1]);
  }

  outputFile << "Finishing up: ------" << std::endl;
  CodecUtils::outputMat(expected_list_rank_1, outputFile);
  outputFile << std::endl;
  CodecUtils::outputMat(expected_list_rank_2, outputFile);
  outputFile << std::endl;
  outputFile << std::endl;

  outputFile.close();
  return 0;
}
