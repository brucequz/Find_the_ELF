#ifndef DUALLISTDECODER_H
#define DUALLISTDECODER_H

#include <map>
#include <queue>
#include <vector>
#include <iostream>
#include <fstream>
#include <climits>
#include <chrono>

#include "minHeap.h"
#include "dualListMap.h"
#include "stopWatch.h"

class FeedForwardTrellis;

struct Cell {
  bool init = false;
  double pathMetric = 3000;
  int fatherState = -1;
  double subPathMetric = 3000;
  int subFatherState = -1;
};

struct MessageInformation {
  MessageInformation(){};
  int decoder_index = -1;
  std::vector<int> message;
  std::vector<int> path;
  std::pair<int, int> begin_end_states;
  double path_metric;
  int list_rank;
  int crc_passing_rank;
  bool list_size_exceeded = false;  // added, used in DLD
};

struct CodeInformation {
  int k;              // input length
  int n;              // output length
  int v;              // memory elements
  int list_size = 1;  // list decoder list size
  int crc_dec;
  int crc_length;
  std::vector<int> generator_poly;
};

namespace dualdecoderutils {

std::vector<int> convertIntToBits(int integer, const int& length);
int hammingDistance(const std::vector<int> x, const std::vector<int>& y);
double euclideanDistance(const std::vector<double>& x,
                         const std::vector<int>& y);
std::vector<int> xOR(const std::vector<int>& x, const std::vector<int>& y);
template <typename T>
void print(const std::vector<T>& vec) {
  for (const T& element : vec) {
    std::cout << element << " ";
  }
  std::cout << std::endl;
}

template <typename T>
void print(const std::vector<std::vector<T>>& matrix) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      std::cout << element << " ";
    }
    std::cout << std::endl;
  }
}

template <typename T>
void output(const std::vector<T>& vec, std::ofstream& outputFile) {
  for (const T& element : vec) {
    outputFile << element << " ";
  }
  // outputFile << std::endl;
}

template <typename T>
void output(const std::vector<std::vector<T>>& matrix,
            std::ofstream& outputFile) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      outputFile << element << " ";
    }
    outputFile << std::endl;
  }
}

template <typename T>
void outputMat(const std::vector<T>& vec, std::ofstream& outputFile) {
  outputFile << "[";
  for (const T& ele : vec) {
    outputFile << ele << ", ";
  }
  outputFile << "] with size: " << vec.size();
}

template <typename T>
bool areVectorsEqual(const std::vector<T>& vector1,
                     const std::vector<T>& vector2) {
  if (vector1.size() != vector2.size()) {
    return false;  // Vectors have different sizes, so they cannot be equal.
  }

  for (size_t i = 0; i < vector1.size(); ++i) {
    if (vector1[i] != vector2[i]) {
      return false;  // Elements at index i are different.
    }
  }

  return true;  // All elements are equal.
}

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
             const std::vector<MessageInformation>& vec2);
}


class DualListDecoder {
 public:
  DualListDecoder(std::vector<CodeInformation> code_info, int max_searched_path);
  DualListDecoder(CodeInformation encoder, std::vector<CodeInformation> code_info, int max_searched_path);
  ~DualListDecoder();
  
  // RATE 1/1 DLD DECODERS

  // Decoding function that alternates between two dual list decoders.
  // This function does not take crc degrees into consideration
  DLDInfo AdaptiveDecode_SimpleAlternate(
      std::vector<double> received_signal, std::vector<std::chrono::microseconds>& timeDurations);
  
  // Decoding function that alternates between two dual list decoders considering crc degrees.
  // For example, if crc degree for two decoders are 3 and 5 respectively. Then every 8 and 32 codewords
  // in both lists there exist one codeword that passes crc (?).
  DLDInfo AdaptiveDecode_CRCAlternate(std::vector<double> received_signal, std::vector<std::chrono::microseconds>& timeDurations);
  
  // RATE 1/2 DLD DECODERS

  // Decoding dunction that contains two rate 1/2 decoders. It alternates between two dual list decoders.
  // This function does not take crc degrees into consideration
  DLDInfo AdaptiveDecode_SimpleAlternate_rate_1_2(std::vector<double> received_signal, std::vector<std::chrono::microseconds>& timeDurations);

  DLDInfo LookAheadDecode_SimpleAlternate_rate_1_2(
      std::vector<double> received_signal, std::vector<std::chrono::microseconds>& timeDurations, std::vector<double> metric_0, std::vector<double> metric_1);
  
  MessageInformation TraceBack_Single(MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
                             const std::vector<std::vector<Cell>>& trellis_states, std::vector<std::vector<int>>& prev_paths,
                             int& num_path_searched, int num_total_stages);
  
  // Trace back according to crc degree
  std::vector<MessageInformation> TraceBack_Multiple(
                            MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
                            const std::vector<std::vector<Cell>>& trellis_states,
                            std::vector<std::vector<int>>& prev_paths, int& num_path_searched,
                            int num_total_stages, int num_trace_back);


  // RATE 1/2 DLD Utilities
  std::vector<double> HardDecode(const std::vector<double>& received_signal);

 private:
  CodeInformation encoder_;
  FeedForwardTrellis* encoder_trellis_ptr_;
  std::vector<CodeInformation> code_info_;
  std::vector<FeedForwardTrellis*> trellis_ptrs_;
  int max_path_to_search_;
  int crc_ratio_; // ratio between crc degree of both dld decoders. For example, crc1 = 3, crc2 = 5. Then crc_ratio = 2^5 / 2^3 = 4;

  //// ZTCC
  // Construct a ZTCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a regular euclidean metric.
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_EuclideanMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::microseconds& ssv_time);

  // Construct a ZTCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a special metric shown by Bill Ryan. Instead of calculating euclidean distance between received point and +/- 1,
  // we compute the product of these two values.
  std::vector<std::vector<Cell>> ConstructZTCCTrellis_WithList_ProductMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::microseconds& ssv_time);

  //// TBCC
  // @todo
  // Construct a TBCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a regular euclidean metric.
  std::vector<std::vector<Cell>> ConstructTBCCTrellis_WithList_EuclideanMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::microseconds& ssv_time);

  // @todo
  // Construct a TBCC Trellis measuring the time taken by trellis construction (for both lists, iteratively)
  // Uses a special metric shown by Bill Ryan. Instead of calculating euclidean distance between received point and +/- 1,
  // we compute the product of these two values.
  std::vector<std::vector<Cell>> ConstructTBCCTrellis_WithList_ProductMetric(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::microseconds& ssv_time);
  
  std::vector<int> convertPathtoMessage(
    const std::vector<int> path, FeedForwardTrellis* trellis_ptr);
  std::vector<int> convertPathtoTrimmedMessage(
    const std::vector<int> path, CodeInformation code, FeedForwardTrellis* trellis_ptr);
  std::vector<int> deconvolveCRC(const std::vector<int>& output, CodeInformation code);
  bool CRC_Check(std::vector<int> input_data, int crc_bits_num, int crc_dec);
};

#endif