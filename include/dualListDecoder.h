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

// struct Cell;
// struct MessageInformation;
// struct CodeInformation;
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

// struct DLDInfo {
//   double combined_metric;
//   std::vector<int> message;
//   std::vector<int> list_ranks;
// };

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

// template <typename T>
// bool areVectorsEqual(const std::vector<T>& vector1, const std::vector<T>& vector2) {
//   if (vector1.size() != vector2.size()) {
//       return false; // Vectors have different sizes, so they cannot be equal.
//   }

//   for (size_t i = 0; i < vector1.size(); ++i) {
//       if (vector1[i] != vector2[i]) {
//           return false; // Elements at index i are different.
//       }
//   }

//   return true; // All elements are equal.
// }

}  // namespace dualdecoderutils



// namespace bpsk {

// std::vector<int> modulate(std::vector<int> encoded_msg) {
//   std::vector<int> modulated_signal(encoded_msg.size());
//   for (int i = 0; i < encoded_msg.size(); ++i) {
//     modulated_signal[i] = -2 * encoded_msg[i] + 1;
//   }
//   return modulated_signal;
// }

// std::vector<int> demodulate(std::vector<double> received_signal) {
//   std::vector<int> hard_decision_demodulated_signal(received_signal.size());
//   for (int i = 0; i < received_signal.size(); ++i) {
//     hard_decision_demodulated_signal[i] = (received_signal[i] < 0.0);
//   }
//   return hard_decision_demodulated_signal;
// }

// }  // namespace bpsk

// namespace crc {

// int binSum(const int& x, const int& y) { return (x + y) % 2; }

// std::vector<int> decToBin(int input, int bit_number) {
//   std::vector<int> output(bit_number, 0);
//   for (int i = bit_number - 1; i >= 0; --i) {
//     output[bit_number - 1 - i] = ((input >> i) & 1);
//   }
//   return output;
// }

// }  // namespace crc

class FeedForwardTrellis;

class DualListDecoder {
 public:
  DualListDecoder(std::vector<CodeInformation> code_info, int max_searched_path);
  ~DualListDecoder();

  DLDInfo adaptiveDecode(
      std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations);
  MessageInformation traceBack(MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
                             const std::vector<std::vector<Cell>>& trellis_states, std::vector<std::vector<int>>& prev_paths,
                             int& num_path_searched, int num_total_stages);

 private:
  std::vector<CodeInformation> code_info_;
  std::vector<FeedForwardTrellis*> trellis_ptrs_;
  int max_path_to_search_;

  std::vector<std::vector<Cell>> constructZTCCTrellis(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr); // deprecated
  std::vector<std::vector<Cell>> constructZTListTrellis(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);
  std::vector<std::vector<Cell>> constructZTListTrellis_precompute(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time);
  MessageInformation ztListDecoding(std::vector<double> receivedMessage);
  std::vector<int> convertPathtoMessage(
    const std::vector<int> path, FeedForwardTrellis* trellis_ptr);
  std::vector<int> convertPathtoTrimmedMessage(
    const std::vector<int> path, CodeInformation code, FeedForwardTrellis* trellis_ptr);
  std::vector<int> deconvolveCRC(const std::vector<int>& output, CodeInformation code);
  bool checkCRC(std::vector<int> demodulated, CodeInformation code);
  bool crc_check(std::vector<int> input_data, int crc_bits_num, int crc_dec);
};

#endif