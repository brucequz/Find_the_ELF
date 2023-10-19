#ifndef DUALLISTDECODER_H
#define DUALLISTDECODER_H

#include <map>
#include <queue>
#include <vector>

#include "minHeap.h"

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

}  // namespace dualdecoderutils

namespace bpsk {

std::vector<int> modulate(std::vector<int> encoded_msg) {
  std::vector<int> modulated_signal(encoded_msg.size());
  for (int i = 0; i < encoded_msg.size(); ++i) {
    modulated_signal[i] = -2 * encoded_msg[i] + 1;
  }
  return modulated_signal;
}

std::vector<int> demodulate(std::vector<double> received_signal) {
  std::vector<int> hard_decision_demodulated_signal(received_signal.size());
  for (int i = 0; i < received_signal.size(); ++i) {
    hard_decision_demodulated_signal[i] = (received_signal[i] < 0.0);
  }
  return hard_decision_demodulated_signal;
}

}  // namespace bpsk

namespace crc {

int binSum(const int& x, const int& y) { return (x + y) % 2; }

std::vector<int> decToBin(int input, int bit_number) {
  std::vector<int> output(bit_number, 0);
  for (int i = bit_number - 1; i >= 0; --i) {
    output[bit_number - 1 - i] = ((input >> i) & 1);
  }
  return output;
}

}  // namespace crc

struct Cell {
  bool init = false;
  double pathMetric = INT_MAX;
  int fatherState = -1;
  double subPathMetric = INT_MAX;
  int subFatherState = -1;
};

struct MessageInformation {
  MessageInformation(){};

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

class FeedForwardTrellis;

class DualListDecoder {
 public:
  DualListDecoder(std::vector<CodeInformation> code_info);
  ~DualListDecoder(){};

  std::vector<MessageInformation> adaptiveDecode(
      std::vector<double> received_signal);
  MessageInformation traceBack(MinHeap* heap, CodeInformation code, FeedForwardTrellis* trellis_ptr,
                             std::vector<std::vector<Cell>> trellis_states, std::vector<std::vector<int>>& prev_paths,
                             int& num_path_searched, int num_total_stages);

 private:
  std::vector<CodeInformation> code_info_;
  std::vector<MinHeap> min_heaps_;
  std::vector<FeedForwardTrellis*> trellis_ptrs_;

  std::vector<std::vector<Cell>> constructZTCCTrellis(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr);
  std::vector<int> convertPathtoMessage(
    const std::vector<int> path, FeedForwardTrellis* trellis_ptr);
  std::vector<int> convertPathtoTrimmedMessage(
    const std::vector<int> path, CodeInformation code, FeedForwardTrellis* trellis_ptr);
  std::vector<int> deconvolveCRC(const std::vector<int>& output, CodeInformation code);
  bool DualListDecoder::checkCRC(std::vector<int> demodulated, CodeInformation code);
};

#endif