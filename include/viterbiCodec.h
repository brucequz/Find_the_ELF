#ifndef VITERBICODEC_H
#define VITERBICODEC_H

#include <climits>
#include <vector>
#include <iostream>
#include <fstream>

class FeedForwardTrellis;

namespace CodecUtils {
std::vector<int> convertIntToBits(int integer, const int& length);
int hammingDistance(const std::vector<int> x, const std::vector<int>& y);
double euclideanDistance(const std::vector<double>& x, const std::vector<int>& y);
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
void output(const std::vector<std::vector<T>>& matrix, std::ofstream& outputFile) {
  for (const std::vector<T>& row : matrix) {
      for (const T& element : row) {
          outputFile << element << " ";
      }
      outputFile << std::endl;
  }
}
}  // namespace CodecUtils

namespace BPSK {
  
std::vector<int> modulate(std::vector<int> encoded_msg);
std::vector<int> demodulate(std::vector<double> received_signal);

} // namespace BPSK

namespace CRC {
  int binSum(const int& x, const int& y);
  std::vector<int> decToBin(int input, int bit_number);
}  // namespace CRC

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
  int k;  // input length
  int n;  // output length
  int v;  // memory elements
  int list_size = 1; // list decoder list size
  int crc_dec;
  int crc_length;
  std::vector<int> generator_poly;
};

class ViterbiCodec {
 public:
  ViterbiCodec(int k, int n, int v, std::vector<int> poly);
  ViterbiCodec(CodeInformation code);
  ~ViterbiCodec();
  std::vector<int> encode(const std::vector<int>& message);
  std::vector<int> encodeZTCC(std::vector<int> message);
  MessageInformation viterbiDecode(const std::vector<int>& coded);
  MessageInformation softViterbiDecode(const std::vector<double>& coded);

  std::vector<int> calculateCRC(const std::vector<int>& input);
  std::vector<int> convolveCRC(const std::vector<int>& input);

  // Function definition in listDecoder.cpp
  std::vector<MessageInformation> listViterbiDecoding(
      const std::vector<double>& received_signal);
  std::vector<MessageInformation> unconstraintListViterbiDecoding(
    const std::vector<double>& received_signal);

 private:
  int k_;  // input message length
  int n_;  // output message length
  int v_;  // number of memory elements
  int crc_dec_; // crc poly in decimal representation
  int crc_length_; // length of crc in binary representation
  double code_rate_;
  int numStates_;
  int list_size_;
  FeedForwardTrellis* trellis_ptr_;

  std::vector<std::vector<Cell>> constructTrellis(
      const std::vector<int>& coded);
  std::vector<std::vector<Cell>> constructTrellis(
      const std::vector<double>& received_signal);
  std::vector<int> convertPathtoMessage(const std::vector<int> path);
  bool checkCRC(std::vector<int> demodulated);
};

#endif