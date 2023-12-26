#include "../include/viterbiCodec.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>

#include "../include/feedForwardTrellis.h"
#include "../include/minHeap.h"
#include "../include/stopWatch.h"

namespace {

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
void print_best_metrics(const std::vector<std::vector<T>>& matrix) {
  for (const std::vector<T>& row : matrix) {
    for (const T& element : row) {
      std::cout << element.pathMetric << " ";
    }
    std::cout << std::endl;
  }
}
}  // namespace
namespace CodecUtils {

std::vector<int> convertIntToBits(int integer, const int& length) {
  if (integer < 0) {
    std::cerr << "CANNOT CONVERT: negative integer" << std::endl;
  } else if (std::ceil(std::log2(integer + 1)) > length) {
    std::cerr << "CANNOT CONVERT: integer too large" << std::endl;
  }
  std::vector<int> result(length, 0);
  int i = length - 1;
  while (integer > 0 && i >= 0) {
    int remainder = integer % 2;
    result[i] = remainder;
    integer /= 2;
    i--;
  }
  return result;
}

int hammingDistance(const std::vector<int> x, const std::vector<int>& y) {
  assert(x.size() == y.size());
  int distance = 0;
  for (int i = 0; i < x.size(); ++i) {
    distance += (x[i] != y[i]);
  }
  return distance;
}

double euclideanDistance(const std::vector<double>& x,
                         const std::vector<int>& y) {
  assert(x.size() == y.size());
  double distance = 0.0;
  for (int i = 0; i < x.size(); ++i) {
    distance += std::pow(x[i] - (double)y[i], 2);
  }
  return distance;
}

std::vector<int> xOR(const std::vector<int>& x, const std::vector<int>& y) {
  std::vector<int> result;
  if (x.size() != y.size()) {
    std::cerr << "INVALID OPERATION SIZE: incompatible operands" << std::endl;
  }
  for (size_t i = 0; i < x.size(); ++i) {
    result.push_back(x[i] ^ y[i]);
  }
  return result;
}

void dec_to_binary(int input, std::vector<int>& output, int bit_number) {
  output.assign(bit_number, -1);
  for (int i = bit_number - 1; i >= 0; i--) {
    int k = input >> i;
    if (k & 1)
      output[bit_number - 1 - i] = 1;
    else
      output[bit_number - 1 - i] = 0;
  }
}

std::vector<int> get_point(int output, int n) {
  std::vector<int> bin_output;
  dec_to_binary(output, bin_output, n);
  for (int i = 0; i < n; i++) {
    bin_output[i] = -2 * bin_output[i] + 1;
  }
  return bin_output;
}

}  // namespace CodecUtils

namespace BPSK {

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

}  // namespace BPSK

namespace CRC {

int binSum(const int& x, const int& y) { return (x + y) % 2; }

std::vector<int> decToBin(int input, int bit_number) {
  std::vector<int> output(bit_number, 0);
  for (int i = bit_number - 1; i >= 0; --i) {
    output[bit_number - 1 - i] = ((input >> i) & 1);
  }
  return output;
}

}  // namespace CRC

ViterbiCodec::ViterbiCodec(int k, int n, int v, std::vector<int> poly)
    : k_(k), n_(n), v_(v) {
  code_rate_ = static_cast<double>(n_ / k_);
  trellis_ptr_ = new FeedForwardTrellis(k, n, v, poly);
  numStates_ = std::pow(2, v);
  list_size_ = 1;
}

ViterbiCodec::ViterbiCodec(CodeInformation code) {
  k_ = code.k;
  n_ = code.n;
  v_ = code.v;
  list_size_ = code.list_size;
  crc_dec_ = code.crc_dec;
  crc_length_ = code.crc_length;
  code_rate_ = static_cast<double>(n_ / k_);
  trellis_ptr_ = new FeedForwardTrellis(k_, n_, v_, code.generator_poly);
  numStates_ = std::pow(2, v_);
}

ViterbiCodec::~ViterbiCodec() {
  delete trellis_ptr_;
  trellis_ptr_ = nullptr;
}

std::vector<int> ViterbiCodec::encode(const std::vector<int>& message) {
  // not necessarily zero terminated

  return trellis_ptr_->encode(message);
}

std::vector<int> ViterbiCodec::encodeZTCC(std::vector<int> message) {
  // append m_ number of zeros to the message
  std::vector<int> encoded_message = trellis_ptr_->encode(message);
  for (int i = 0; i < v_; ++i) {
    message.push_back(0);
  }
  return trellis_ptr_->encode(message);
}

MessageInformation ViterbiCodec::softViterbiDecoding(
    std::vector<double> receivedMessage, std::chrono::milliseconds& ssv_traceback_time) {
  std::vector<std::vector<Cell>> trellisInfo;
  int lowrate_pathLength = (receivedMessage.size() / n_) + 1;

  trellisInfo = std::vector<std::vector<Cell>>(
      numStates_, std::vector<Cell>(lowrate_pathLength));

  // initializes all the valid starting states

  trellisInfo[0][0].pathMetric = 0;
  trellisInfo[0][0].init = true;

  // Add stopwatch
  Stopwatch soft_sw;
  soft_sw.tic();

  // building the trellis
  for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
    for (int currentState = 0; currentState < numStates_; currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init) continue;

      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0;
           forwardPathIndex < trellis_ptr_->nextStates_[0].size();
           forwardPathIndex++) {
        // since our transitions correspond to symbols, the forwardPathIndex has
        // no correlation beyond indexing the forward path

        int nextState =
            trellis_ptr_->nextStates_[currentState][forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0) continue;

        double branchMetric = 0;
        std::vector<int> output_point = CodecUtils::get_point(
            trellis_ptr_->output_[currentState][forwardPathIndex], n_);

        for (int i = 0; i < n_; i++) {
          branchMetric += std::pow(
              receivedMessage[n_ * stage + i] - (double)output_point[i], 2);
          // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
          // stage + i] - (double)output_point[i]);
        }
        double totalPathMetric =
            branchMetric + trellisInfo[currentState][stage].pathMetric;

        // dealing with cases of uninitialized states, when the transition
        // becomes the optimal father state, and suboptimal father state, in
        // order
        if (!trellisInfo[nextState][stage + 1].init) {
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
          trellisInfo[nextState][stage + 1].init = true;
        } else if (trellisInfo[nextState][stage + 1].pathMetric >
                   totalPathMetric) {
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
        }
      }
    }
  }

  // perform the single traceback
  MessageInformation output;

  double minMetric = INT_MAX;
  int endingState = -1;

  if (trellisInfo[0][lowrate_pathLength - 1].pathMetric < minMetric) {
    minMetric = trellisInfo[0][lowrate_pathLength - 1].pathMetric;
    endingState = 0;
  }

  std::vector<int> path(lowrate_pathLength);
  path[lowrate_pathLength - 1] = endingState;
  int currentState = endingState;

  // traceback
  for (int stage = lowrate_pathLength - 1; stage > 0; stage--) {
    currentState = trellisInfo[currentState][stage].fatherState;
    path[stage - 1] = currentState;
  }

  soft_sw.toc();
  ssv_traceback_time += soft_sw.getElapsed();
  soft_sw.reset();
  
  std::vector<int> message = convertPathtoTrimmedMessage(path);
  output.message = message;
  output.path = path;
  output.path_metric = minMetric;
  output.begin_end_states = {path.front(), path.back()};
  return output;
}

std::vector<std::vector<Cell>> ViterbiCodec::constructZTTrellis(
    std::vector<double> receivedMessage) {
  std::vector<std::vector<Cell>> trellisInfo;
  int lowrate_pathLength = (receivedMessage.size() / n_) + 1;
  trellisInfo = std::vector<std::vector<Cell>>(
      numStates_, std::vector<Cell>(lowrate_pathLength));

  // initializes all the valid starting states

  trellisInfo[0][0].pathMetric = 0;
  trellisInfo[0][0].init = true;

  // building the trellis
  for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
    for (int currentState = 0; currentState < numStates_; currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init) continue;

      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0;
           forwardPathIndex < trellis_ptr_->nextStates_[0].size();
           forwardPathIndex++) {
        // since our transitions correspond to symbols, the forwardPathIndex has
        // no correlation beyond indexing the forward path

        int nextState =
            trellis_ptr_->nextStates_[currentState][forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0) continue;

        double branchMetric = 0;
        std::vector<int> output_point = CodecUtils::get_point(
            trellis_ptr_->output_[currentState][forwardPathIndex], n_);

        for (int i = 0; i < n_; i++) {
          branchMetric += std::pow(
              receivedMessage[n_ * stage + i] - (double)output_point[i], 2);
          // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
          // stage + i] - (double)output_point[i]);
        }
        double totalPathMetric =
            branchMetric + trellisInfo[currentState][stage].pathMetric;

        // dealing with cases of uninitialized states, when the transition
        // becomes the optimal father state, and suboptimal father state, in
        // order
        if (!trellisInfo[nextState][stage + 1].init) {
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
          trellisInfo[nextState][stage + 1].init = true;
        } else if (trellisInfo[nextState][stage + 1].pathMetric >
                   totalPathMetric) {
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
        }
      }
    }
  }
  return trellisInfo;
}

std::vector<int> ViterbiCodec::calculateCRC(const std::vector<int>& input) {
  // generating (crc_length - 1) number of redundancy bits (crc bits)
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);

  std::vector<int> output = input;
  output.resize(input.size() + crc_length_ - 1, 0);

  // long division
  for (int i = 0; i < input.size(); ++i) {
    if (output[i] == 1) {
      std::transform(output.begin() + i, output.begin() + i + crc_length_,
                     crc_bin.begin(), output.begin() + i, CRC::binSum);
    }
  }

  std::copy(input.begin(), input.end(), output.begin());

  return output;
}

std::vector<int> ViterbiCodec::convolveCRC(const std::vector<int>& input) {
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);
  int input_size = input.size();
  int crc_size = crc_bin.size();
  int output_size = input_size + crc_size - 1;
  std::vector<int> output(output_size, 0);

  for (int i = 0; i < input_size; ++i) {
    for (int j = 0; j < crc_size; ++j) {
      output[i + j] += input[i] * crc_bin[j];
    }
  }

  std::for_each(output.begin(), output.end(),
                [](int& element) { element %= 2; });

  return output;
}

std::vector<int> ViterbiCodec::deconvolveCRC(const std::vector<int>& output) {
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);
  std::vector<int> result(output.size() - crc_bin.size() + 1, 0);
  for (uint64_t n = 0; n < result.size(); n++) {
    result[n] = output[n];
    uint64_t start = std::max((int)(n - crc_bin.size() + 1), 0);
    for (uint64_t i = start; i < n; i++) {
      result[n] -= result[i] * crc_bin[n - i];
    }
    result[n] = (result[n] >= 0) ? result[n] : -result[n];
    result[n] /= crc_bin[0];
  }
  std::for_each(result.begin(), result.end(),
                [](int& element) { element %= 2; });
  return result;
}

bool ViterbiCodec::checkCRC(std::vector<int> demodulated) {
  // check crc by dividing the demodulated signal with crc poly
  std::vector<int> crc_bin = CRC::decToBin(crc_dec_, crc_length_);

  for (int ii = 0; ii <= (int)demodulated.size() - crc_length_; ii++) {
    if (demodulated[ii] == 1) {
      // Note: transform doesn't include .end
      std::transform(demodulated.begin() + ii,
                     demodulated.begin() + (ii + crc_length_), crc_bin.begin(),
                     demodulated.begin() + ii, CRC::binSum);
    }
  }
  bool all_zero = std::all_of(demodulated.begin(), demodulated.end(),
                              [](int i) { return i == 0; });
  return all_zero;
}