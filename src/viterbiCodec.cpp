#include "../include/viterbiCodec.h"

#include <cassert>
#include <cmath>
#include <iostream>

#include "../include/feedForwardTrellis.h"
#include "../include/minHeap.h"

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
  for (int i = 0; i < v_; ++i) {
    message.push_back(0);
  }
  return trellis_ptr_->encode(message);
}

MessageInformation ViterbiCodec::viterbiDecode(const std::vector<int>& coded) {
  std::vector<std::vector<Cell>> trellis_states = constructTrellis(coded);
  MessageInformation output;

  int num_total_stages = trellis_states[0].size();
  MinHeap heap;  // Detour Tree

  // add all final stage nodes to the heap
  for (int i = 0; i < numStates_; ++i) {
    DetourNode node;
    node.start_state = i;
    node.path_metric = trellis_states[i][num_total_stages - 1].pathMetric;
    heap.insert(node);
  }

  // pop the node with the smallest path metric - ML
  DetourNode min_node = heap.pop();
  std::vector<int> path(num_total_stages);
  std::vector<int> message(k_ * (num_total_stages - 1), 0);

  int cur_state = min_node.start_state;

  for (int stage = num_total_stages - 1; stage >= 0; --stage) {
    int father_state = trellis_states[cur_state][stage].fatherState;
    path[stage] = cur_state;
    if (stage == 0) {
      break;
    }
    // assuming k_ == 1
    assert(k_ == 1);
    int input =
        (cur_state == trellis_ptr_->nextStates_[father_state][0]) ? 0 : 1;
    message[stage - 1] = input;
    cur_state = father_state;
  }
  output.message = message;
  output.path = path;
  return output;
}

MessageInformation ViterbiCodec::softViterbiDecode(
    const std::vector<double>& received_signal) {
  MessageInformation output;
  std::vector<std::vector<Cell>> trellis_states =
      constructTrellis(received_signal);

  int num_total_stages = trellis_states[0].size();
  MinHeap heap;  // Detour Tree

  // add all final stage nodes to the heap
  for (int i = 0; i < numStates_; ++i) {
    DetourNode node;
    node.start_state = i;
    node.path_metric = trellis_states[i][num_total_stages - 1].pathMetric;
    heap.insert(node);
  }

  // pop the node with the smallest path metric - ML
  DetourNode min_node = heap.pop();
  if (min_node.start_state != 0) {
    min_node = heap.pop();
  }
  std::vector<int> path(num_total_stages);
  std::vector<int> message(k_ * (num_total_stages - 1), 0);

  std::cout << "Min path metric: " << min_node.path_metric << std::endl;

  int cur_state = min_node.start_state;

  for (int stage = num_total_stages - 1; stage >= 0; --stage) {
    int father_state = trellis_states[cur_state][stage].fatherState;
    path[stage] = cur_state;
    if (stage == 0) {
      break;
    }
    // assuming k_ == 1
    assert(k_ == 1);
    int input =
        (cur_state == trellis_ptr_->nextStates_[father_state][0]) ? 0 : 1;
    message[stage - 1] = input;
    cur_state = father_state;
  }
  output.message = message;
  output.path = path;
  return output;
}

std::vector<std::vector<Cell>> ViterbiCodec::constructTrellis(
    const std::vector<int>& coded) {
  /*
   *  Construct a reverse trellis states with cells containing relevant
   * information for viterbi decoding.
   *  Cell
   * {
   *    - branchMetric:
   *    - pathMetric:
   *    - init: whether or not a cell is activated, aka. reachable in the
   * decoding process.
   * }
   */

  std::vector<std::vector<Cell>> trellis_states;
  int message_length = coded.size();
  int total_length = message_length / n_;
  trellis_states.resize(numStates_, std::vector<Cell>(total_length + 1));
  trellis_states[0][0].pathMetric = 0;
  trellis_states[0][0].init = true;

  for (int cur_stage = 0; cur_stage < total_length; ++cur_stage) {
    for (int cur_state = 0; cur_state < numStates_; ++cur_state) {
      if (!trellis_states[cur_state][cur_stage].init) {
        continue;
      }
      int cur_path_metric = trellis_states[cur_state][cur_stage].pathMetric;
      auto begin = coded.begin() + cur_stage * n_;
      auto end = begin + n_;
      std::vector<int> target_message(begin, end);
      // activate the next states
      for (int i = 0; i < trellis_ptr_->nextStates_[cur_state].size(); ++i) {
        int next_state = trellis_ptr_->nextStates_[cur_state][i];
        // trellis_states[next_state][cur_stage + 1].init = true;

        int possible_output = trellis_ptr_->output_[cur_state][i];
        std::vector<int> expected_output =
            CodecUtils::convertIntToBits(possible_output, n_);

        int branch_metric =
            CodecUtils::hammingDistance(expected_output, target_message);
        int temp_path_metric = cur_path_metric + branch_metric;

        Cell* target_cell = &trellis_states[next_state][cur_stage + 1];

        if (!target_cell->init) {
          // if the next state is not initialized, we temporarily store the path
          // metric
          target_cell->init = true;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else if (target_cell->pathMetric > temp_path_metric) {
          // the current path metric is better
          target_cell->subPathMetric = target_cell->pathMetric;
          target_cell->subFatherState = target_cell->fatherState;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else {
          // the current path metric is worse
          target_cell->subPathMetric = temp_path_metric;
          target_cell->subFatherState = cur_state;
        }
      }
    }
  }
  return trellis_states;
}

std::vector<std::vector<Cell>> ViterbiCodec::constructTrellis(
    const std::vector<double>& received_signal) {
  std::vector<std::vector<Cell>> trellis_states;
  int message_length = received_signal.size();
  int total_length = message_length / n_;
  trellis_states.resize(numStates_, std::vector<Cell>(total_length + 1));
  trellis_states[0][0].pathMetric = 0.0;
  trellis_states[0][0].init = true;

  for (int cur_stage = 0; cur_stage < total_length; ++cur_stage) {
    for (int cur_state = 0; cur_state < numStates_; ++cur_state) {
      if (!trellis_states[cur_state][cur_stage].init) {
        continue;
      }
      double cur_path_metric = trellis_states[cur_state][cur_stage].pathMetric;
      auto begin = received_signal.begin() + cur_stage * n_;
      auto end = begin + n_;
      std::vector<double> target_message(begin, end);

      // activate the next states
      for (int i = 0; i < trellis_ptr_->nextStates_[cur_state].size(); ++i) {
        int next_state = trellis_ptr_->nextStates_[cur_state][i];
        // trellis_states[next_state][cur_stage + 1].init = true;

        int possible_output = trellis_ptr_->output_[cur_state][i];
        std::vector<int> expected_output =
            CodecUtils::convertIntToBits(possible_output, n_);

        // modualted expected output
        std::vector<int> expected_signal = BPSK::modulate(expected_output);

        double branch_metric =
            CodecUtils::euclideanDistance(target_message, expected_signal);
        double temp_path_metric = cur_path_metric + branch_metric;

        Cell* target_cell = &trellis_states[next_state][cur_stage + 1];

        if (!target_cell->init) {
          // if the next state is not initialized, we temporarily store the path
          // metric
          target_cell->init = true;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else if (target_cell->pathMetric > temp_path_metric) {
          // the current path metric is better
          target_cell->subPathMetric = target_cell->pathMetric;
          target_cell->subFatherState = target_cell->fatherState;
          target_cell->pathMetric = temp_path_metric;
          target_cell->fatherState = cur_state;
        } else {
          // the current path metric is worse
          target_cell->subPathMetric = temp_path_metric;
          target_cell->subFatherState = cur_state;
        }
      }
    }
  }
  return trellis_states;
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