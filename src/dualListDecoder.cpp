#include "../include/dualListDecoder.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <vector>

#include "../include/dualListMap.h"
#include "../include/feedForwardTrellis.h"

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
}  // namespace
namespace dualdecoderutils {
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

DualListDecoder::DualListDecoder(std::vector<CodeInformation> code_info,
                                 int max_searched_path)
    : code_info_(code_info), max_path_to_search_(max_searched_path) {
  CodeInformation code_list_0 = code_info[0];
  CodeInformation code_list_1 = code_info[1];

  FeedForwardTrellis* trellis_ptr_0 = new FeedForwardTrellis(code_list_0);
  FeedForwardTrellis* trellis_ptr_1 = new FeedForwardTrellis(code_list_1);
  trellis_ptrs_.push_back(trellis_ptr_0);
  trellis_ptrs_.push_back(trellis_ptr_1);
}

DualListDecoder::~DualListDecoder() {
  for (FeedForwardTrellis*& trellis_ptr : trellis_ptrs_) {
    delete trellis_ptr;
    trellis_ptr = nullptr;
  }
}

DLDInfo DualListDecoder::adaptiveDecode(std::vector<double> received_signal) {
  /*
  This function adaptively expand the list size until the smallest future match
  (SFM) metric is larger than the best current match (BCM) metric.

  local variables:
    - double best_current_match;
    - double smallest_future_match;
       - smallest metric on the left +

  caveat: cache miss (due to size of minHeap)

  algorithm adaptive decoder:
    Input: a received signal and code polynomials and crc polynomials
    Output: a vector of agreed messages (maximum likelihood message)

      create an unordered map to store available messages
      while ML path is not found
        add a new l1 path to the map
          if path already exists in map
            end while
          end if
        add a new l2 path to the map
          if path already exists in map
            end while
          end if
      end while

      if the most likely agreed message is found
        set BCM to the combined path metric
        create thresholds constraint on both list decoders. (l1, l2)

        if (a new path generated on l1 >= l1 threshold)
          stop path tracing from l1
        end if
        if (a new path generated on l2 >= l2 threshold)
          stop path tracing from l2
        end if

        if both list decoders are stopped
          return maximum likelihood message
        end if
      end if
  */
  std::vector<std::vector<MessageInformation>> output;
  std::vector<MessageInformation> output_0;
  std::vector<MessageInformation> output_1;
  DualListMap mp;
  double best_current_match = INT_MAX;
  double decoder_threshold_0 = INT_MAX;
  double decoder_threshold_1 = INT_MAX;

  // divide the received signal
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

  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      constructZTCCTrellis(received_codec_1, code_0, trellis_ptrs_[0]);
  int num_total_stages_0 = trellis_0[0].size();
  std::vector<std::vector<int>> prev_paths_list_0;
  MinHeap* heap_list_0 = new MinHeap;

  DetourNode node_0;
  node_0.start_state = 0;
  node_0.path_metric = trellis_0[0][num_total_stages_0 - 1].pathMetric;
  heap_list_0->insert(node_0);
  int num_path_searched_0 = 0;
  bool decoder_0_stop = false;

  // list decoder 1
  std::vector<std::vector<Cell>> trellis_1 =
      constructZTCCTrellis(received_codec_2, code_1, trellis_ptrs_[1]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;

  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 > max_path_to_search_) {
      decoder_0_stop = true;
    }
    if (!decoder_0_stop) {
      // std::cout << "decoder 0 traceback" << std::endl;
      MessageInformation mi_0 =
          traceBack(heap_list_0, code_0, trellis_ptrs_[0], trellis_0,
                    prev_paths_list_0, num_path_searched_0, num_total_stages_0);
      mi_0.decoder_index = 0;
      if (mi_0.path_metric >= decoder_threshold_0) {
        decoder_0_stop = true;
      }
      mp.insert(mi_0);
      output_0.push_back(mi_0);
    }
    // list decoder 1 traceback
    if (num_path_searched_1 > max_path_to_search_) {
      decoder_1_stop = true;
    }
    if (!decoder_1_stop) {
      // std::cout << "decoder 1 traceback" << std::endl;
      MessageInformation mi_1 =
          traceBack(heap_list_1, code_1, trellis_ptrs_[1], trellis_1,
                    prev_paths_list_1, num_path_searched_1, num_total_stages_1);
      mi_1.decoder_index = 1;
      if (mi_1.path_metric >= decoder_threshold_1) {
        decoder_1_stop = true;
      }
      mp.insert(mi_1);
      output_1.push_back(mi_1);
    }

    if (decoder_0_stop && decoder_1_stop) {
      break;
    }

    if (mp.queue_size() != 0) {
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      decoder_threshold_0 = best_current_match - node_1.path_metric;
      decoder_threshold_1 = best_current_match - node_0.path_metric;
      best_combined_found = true;
      // free pointers
      delete heap_list_0;
      delete heap_list_1;
      heap_list_0 = nullptr;
      heap_list_1 = nullptr;
      // std::cout << "found agreed message" << std::endl;
      return agreed_message;
    }
    // std::cout << "hello" << std::endl;
  }

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  empty_message.message = std::vector<int>(64, -1);

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  // std::cout << "found nothing!" << std::endl;
  return empty_message;
}

MessageInformation DualListDecoder::traceBack(
    MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
    const std::vector<std::vector<Cell>>& trellis_states,
    std::vector<std::vector<int>>& prev_paths, int& num_path_searched,
    int num_total_stages) {
  MessageInformation mi;
  bool found_path = false;
  while (!found_path) {
    DetourNode detour = heap->pop();
    std::vector<int> path(num_total_stages);

    int resume_stage = num_total_stages - 1;
    double forward_partial_path_metric = 0.0;
    int cur_state = detour.start_state;

    if (detour.original_path != -1) {
      forward_partial_path_metric = detour.forward_path_metric;
      resume_stage = detour.detour_stage;

      path = prev_paths[detour.original_path];
      cur_state = path[resume_stage];

      double cur_sub_path_metric =
          trellis_states[cur_state][resume_stage].subPathMetric;

      cur_state = trellis_states[cur_state][resume_stage].subFatherState;

      resume_stage--;
      double prev_path_metric =
          trellis_states[cur_state][resume_stage].pathMetric;
      forward_partial_path_metric += cur_sub_path_metric - prev_path_metric;
    }

    path[resume_stage] = cur_state;

    for (int stage = resume_stage; stage > 0; stage--) {
      double cur_sub_path_metric =
          trellis_states[cur_state][stage].subPathMetric;
      double cur_path_metric = trellis_states[cur_state][stage].pathMetric;

      if (trellis_states[cur_state][stage].subFatherState != -1) {
        DetourNode localDetour;
        localDetour.start_state = 0;
        localDetour.path_metric =
            cur_sub_path_metric + forward_partial_path_metric;
        localDetour.forward_path_metric = forward_partial_path_metric;
        localDetour.original_path = num_path_searched;
        localDetour.detour_stage = stage;
        heap->insert(localDetour);
      }

      cur_state = trellis_states[cur_state][stage].fatherState;
      double prev_path_metric = trellis_states[cur_state][stage - 1].pathMetric;
      forward_partial_path_metric += cur_path_metric - prev_path_metric;
      path[stage - 1] = cur_state;
    }
    prev_paths.push_back(path);

    std::vector<int> full_message = convertPathtoMessage(path, trellis_ptr);
    std::vector<int> messageWithoutTrailingZeros =
        convertPathtoTrimmedMessage(path, code, trellis_ptr);
    std::vector<int> message = deconvolveCRC(messageWithoutTrailingZeros, code);

    if (checkCRC(messageWithoutTrailingZeros, code) &&
        path.front() == path.back() && path.back() == 0) {
      mi.path = path;
      mi.path_metric = detour.path_metric;
      // convert path to message
      mi.message = message;
      mi.list_rank = num_path_searched;
      found_path = true;
    }

    num_path_searched++;
  }
  return mi;
}

std::vector<std::vector<Cell>> DualListDecoder::constructZTCCTrellis(
    const std::vector<double>& received_signal, CodeInformation code,
    FeedForwardTrellis* trellis_ptr) {
  /*
  Construct output trellis for traceback later with euclidean distances
  as path metrics.
  */

  std::vector<std::vector<Cell>> trellis_states;
  int signal_length = received_signal.size();
  int number_of_stages = signal_length / code.n;

  trellis_states.resize(std::pow(2, code.v),
                        std::vector<Cell>(number_of_stages + 1));
  trellis_states[0][0].init = true;
  trellis_states[0][0].pathMetric = 0.0;

  for (int cur_stage = 0; cur_stage < number_of_stages; ++cur_stage) {
    for (int cur_state = 0; cur_state < std::pow(2, code.v); ++cur_state) {
      if (!trellis_states[cur_state][cur_stage].init) {
        continue;
      }
      double cur_path_metric = trellis_states[cur_state][cur_stage].pathMetric;
      auto begin = received_signal.begin() + cur_stage * code.n;
      auto end = begin + code.n;
      std::vector<double> target_message(begin, end);

      // activate the next states
      for (int i = 0; i < trellis_ptr->nextStates_[cur_state].size(); ++i) {
        int next_state = trellis_ptr->nextStates_[cur_state][i];
        // trellis_states[next_state][cur_stage + 1].init = true;

        int possible_output = trellis_ptr->output_[cur_state][i];
        std::vector<int> expected_output =
            dualdecoderutils::convertIntToBits(possible_output, code.n);

        // modualted expected output
        std::vector<int> expected_signal = BPSK::modulate(expected_output);

        double branch_metric = dualdecoderutils::euclideanDistance(
            target_message, expected_signal);
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

std::vector<int> DualListDecoder::convertPathtoMessage(
    const std::vector<int> path, FeedForwardTrellis* trellis_ptr) {
  std::vector<int> message;
  for (int path_id = 0; path_id < path.size() - 1; path_id++) {
    for (int forward_path = 0;
         forward_path < trellis_ptr->nextStates_[0].size(); forward_path++) {
      if (trellis_ptr->nextStates_[path[path_id]][forward_path] ==
          path[path_id + 1]) {
        message.push_back(forward_path);
      }
    }
  }
  return message;
}

std::vector<int> DualListDecoder::convertPathtoTrimmedMessage(
    const std::vector<int> path, CodeInformation code,
    FeedForwardTrellis* trellis_ptr) {
  std::vector<int> message;
  message = convertPathtoMessage(path, trellis_ptr);
  message.resize(message.size() - code.v);
  return message;
}

std::vector<int> DualListDecoder::deconvolveCRC(const std::vector<int>& output,
                                                CodeInformation code) {
  std::vector<int> crc_bin = CRC::decToBin(code.crc_dec, code.crc_length);
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

bool DualListDecoder::checkCRC(std::vector<int> demodulated,
                               CodeInformation code) {
  // check crc by dividing the demodulated signal with crc poly
  std::vector<int> crc_bin = CRC::decToBin(code.crc_dec, code.crc_length);

  for (int ii = 0; ii <= (int)demodulated.size() - code.crc_length; ii++) {
    if (demodulated[ii] == 1) {
      // Note: transform doesn't include .end
      std::transform(demodulated.begin() + ii,
                     demodulated.begin() + (ii + code.crc_length),
                     crc_bin.begin(), demodulated.begin() + ii, CRC::binSum);
    }
  }
  bool all_zero = std::all_of(demodulated.begin(), demodulated.end(),
                              [](int i) { return i == 0; });
  return all_zero;
}