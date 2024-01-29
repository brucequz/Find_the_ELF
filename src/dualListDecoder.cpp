/**
 * @file dualListDecoder.cpp
 * @author Bruce Qu (brucequ@ucla.edu)
 * @brief
 * @version 0.1
 * @date 2023-12-27
 *
 * @copyright Copyright (c) 2023
 *
 */
#include "../include/dualListDecoder.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <map>
#include <vector>

#include "../include/dualListMap.h"
#include "../include/feedForwardTrellis.h"
#include "../include/stopWatch.h"

/**
 * @namespace private namespace
 *
 */
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

std::vector<double> ComputeSquaredDifferences(
    const std::vector<double>& vector1, const std::vector<double>& vector2) {
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
    // You can handle this error in your preferred way, e.g., throw an exception
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

/**
 * @namespace dualdecoderutils namespace
 *
 */
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

int bin_sum(int i, int j) { return (i + j) % 2; }
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

  if (code_list_0.v == 0 || code_list_1.v == 0) {
    std::cerr << "Cannot calculate the ratio when one of the numbers is zero."
              << std::endl;
  }
  int largerCRCDegree =
      (code_list_0.v >= code_list_1.v) ? code_list_0.v : code_list_1.v;
  int smallerCRCDegree =
      (code_list_0.v <= code_list_1.v) ? code_list_0.v : code_list_1.v;

  // Calculate the ratio
  crc_ratio_ =
      int(std::pow(2, largerCRCDegree) / std::pow(2, smallerCRCDegree));
  // std::cout << "crc_ratio_" << crc_ratio_ << std::endl;
}

DualListDecoder::DualListDecoder(CodeInformation encoder,
                                 std::vector<CodeInformation> code_info,
                                 int max_searched_path)
    : encoder_(encoder),
      code_info_(code_info),
      max_path_to_search_(max_searched_path) {
  CodeInformation code_list_0 = code_info[0];
  CodeInformation code_list_1 = code_info[1];

  encoder_trellis_ptr_ = new FeedForwardTrellis(encoder_);
  FeedForwardTrellis* trellis_ptr_0 = new FeedForwardTrellis(code_list_0);
  FeedForwardTrellis* trellis_ptr_1 = new FeedForwardTrellis(code_list_1);
  trellis_ptrs_.push_back(trellis_ptr_0);
  trellis_ptrs_.push_back(trellis_ptr_1);

  if (code_list_0.v == 0 || code_list_1.v == 0) {
    std::cerr << "Cannot calculate the ratio when one of the numbers is zero."
              << std::endl;
  }
  int largerCRCDegree =
      (code_list_0.v >= code_list_1.v) ? code_list_0.v : code_list_1.v;
  int smallerCRCDegree =
      (code_list_0.v <= code_list_1.v) ? code_list_0.v : code_list_1.v;

  // Calculate the ratio
  crc_ratio_ =
      int(std::pow(2, largerCRCDegree) / std::pow(2, smallerCRCDegree));
  // std::cout << "crc_ratio_" << crc_ratio_ << std::endl;
}

DualListDecoder::~DualListDecoder() {
  for (FeedForwardTrellis*& trellis_ptr : trellis_ptrs_) {
    delete trellis_ptr;
    trellis_ptr = nullptr;
  }
}

DLDInfo DualListDecoder::adaptiveDecode(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations) {
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
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      constructZTListTrellis_precompute(
          received_codec_1, code_0, trellis_ptrs_[0], timeDurations[0]);
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
      constructZTListTrellis_precompute(
          received_codec_2, code_1, trellis_ptrs_[1], timeDurations[0]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;
  
  // step 1 end: SSV add-compare-select and initial traceback time for both decoders
  // auto ssv_end_time = std::chrono::steady_clock::now();
  // auto ssv_duration = std::chrono::duration_cast<std::chrono::milliseconds>(ssv_end_time - ssv_start_time);

  // step 2: time taken for additional traceback opearations required by SLVD for both decoders
  

  // time taken to do additional traceback and insertion
  Stopwatch sw_step2_3;
  sw_step2_3.tic();
  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 >= max_path_to_search_) {
      decoder_0_stop = true;
    }
    if (!decoder_0_stop) {
      
      // std::cout << "decoder 0 traceback" << std::endl;
      MessageInformation mi_0 =
          TraceBack_Single(heap_list_0, code_0, trellis_ptrs_[0], trellis_0,
                    prev_paths_list_0, num_path_searched_0, num_total_stages_0);
      mi_0.decoder_index = 0;
      if (mi_0.path_metric != -1.0) {
        // if list size is not exceeded
        if (mi_0.path_metric >= decoder_threshold_0) {
          decoder_0_stop = true;
        }
        mp.insert(mi_0);
        output_0.push_back(mi_0);
      }
    }
    // list decoder 1 traceback
    if (num_path_searched_1 >= max_path_to_search_) {
      decoder_1_stop = true;
    }
    if (!decoder_1_stop) {
      // std::cout << "decoder 1 traceback" << std::endl;
      // std::cout << "before " << num_path_searched_1 << std::endl;
      MessageInformation mi_1 =
          TraceBack_Single(heap_list_1, code_1, trellis_ptrs_[1], trellis_1,
                    prev_paths_list_1, num_path_searched_1, num_total_stages_1);
      // std::cout << "after " << num_path_searched_1 << std::endl;
      mi_1.decoder_index = 1;
      if (mi_1.path_metric != -1.0) {
        if (mi_1.path_metric >= decoder_threshold_1) {
          decoder_1_stop = true;
        }
        mp.insert(mi_1);
        output_1.push_back(mi_1);
      }
    }

    if (decoder_0_stop && decoder_1_stop) {
      // std::cout << "Both decoders declared stop" << std::endl;
      break;
    }

    if (mp.queue_size() != 0) {
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // record the received signal just in case the decoding is incorrect
      agreed_message.received_signal = received_signal;
      decoder_threshold_0 = best_current_match - node_1.path_metric;
      decoder_threshold_1 = best_current_match - node_0.path_metric;
      best_combined_found = true;
      // free pointers
      delete heap_list_0;
      delete heap_list_1;
      heap_list_0 = nullptr;
      heap_list_1 = nullptr;
      // std::cout << "found agreed message" << std::endl;
      // std::cout << "Debug: agreed message list size: ";
      // dualdecoderutils::print(agreed_message.list_ranks);
      // std::cout << std::endl;

      return agreed_message;
    }
    // std::cout << "hello" << std::endl;
  }
  sw_step2_3.toc();
  timeDurations[2] += sw_step2_3.getElapsed();
  sw_step2_3.reset();

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  // when the output is not found, we store the message to be all -1's
  // and the received_signal
  empty_message.message = std::vector<int>(64, -1);
  empty_message.received_signal = received_signal;

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  // std::cout << "found nothing!" << std::endl;
  return empty_message;
}

DLDInfo DualListDecoder::AdaptiveDecode_SimpleAlternate(
    std::vector<double> received_signal,
    std::vector<std::chrono::milliseconds>& timeDurations) {
  /**
  @brief This function adaptively expand the list size until the smallest future
  match (SFM) metric is larger than the best current match (BCM) metric.

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
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      constructZTListTrellis_precompute(
          received_codec_1, code_0, trellis_ptrs_[0], timeDurations[0]);
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
      constructZTListTrellis_precompute(
          received_codec_2, code_1, trellis_ptrs_[1], timeDurations[0]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;

  // time taken to do additional traceback and insertion
  Stopwatch sw_step2_3;
  sw_step2_3.tic();
  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 >= max_path_to_search_) {
      decoder_0_stop = true;
    }
    if (!decoder_0_stop) {
      // std::cout << "decoder 0 traceback" << std::endl;
      MessageInformation mi_0 = TraceBack_Single(
          heap_list_0, code_0, trellis_ptrs_[0], trellis_0, prev_paths_list_0,
          num_path_searched_0, num_total_stages_0);
      mi_0.decoder_index = 0;
      if (!mi_0.list_size_exceeded) {
        // if list size is not exceeded
        if (mi_0.path_metric >= decoder_threshold_0) {
          decoder_0_stop = true;
        }
        mp.insert(mi_0);
        output_0.push_back(mi_0);
      }
    }
    // list decoder 1 traceback
    if (num_path_searched_1 >= max_path_to_search_) {
      decoder_1_stop = true;
    }
    if (!decoder_1_stop) {
      // std::cout << "decoder 1 traceback" << std::endl;
      // std::cout << "before " << num_path_searched_1 << std::endl;
      MessageInformation mi_1 = TraceBack_Single(
          heap_list_1, code_1, trellis_ptrs_[1], trellis_1, prev_paths_list_1,
          num_path_searched_1, num_total_stages_1);
      // std::cout << "after " << num_path_searched_1 << std::endl;
      mi_1.decoder_index = 1;
      if (!mi_1.list_size_exceeded) {
        if (mi_1.path_metric >= decoder_threshold_1) {
          decoder_1_stop = true;
        }
        mp.insert(mi_1);
        output_1.push_back(mi_1);
      }
    }

    if (decoder_0_stop && decoder_1_stop) {
      // std::cout << "Both decoders declared stop" << std::endl;
      break;
    }

    if (mp.queue_size() != 0) {
      // std::cout << "Agreed message found: " << std::endl;
      // std::cout << "Simple alternate algorithm ended with path searched: ["
      // << num_path_searched_0 << ", " << num_path_searched_1 << std::endl;
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // record the received signal just in case the decoding is incorrect
      agreed_message.received_signal = received_signal;
      decoder_threshold_0 = best_current_match - node_1.path_metric;
      decoder_threshold_1 = best_current_match - node_0.path_metric;
      best_combined_found = true;
      // free pointers
      delete heap_list_0;
      delete heap_list_1;
      heap_list_0 = nullptr;
      heap_list_1 = nullptr;
      // std::cout << "found agreed message" << std::endl;
      // std::cout << "Debug: agreed message list size: ";
      // dualdecoderutils::print(agreed_message.list_ranks);
      // std::cout << std::endl;

      return agreed_message;
    }
    // std::cout << "hello" << std::endl;
  }
  sw_step2_3.toc();
  timeDurations[2] += sw_step2_3.getElapsed();
  sw_step2_3.reset();

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  // when the output is not found, we store the message to be all -1's
  // and the received_signal
  empty_message.message = std::vector<int>(64, -1);
  empty_message.received_signal = received_signal;

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  // std::cout << "found nothing!" << std::endl;
  return empty_message;
}

DLDInfo DualListDecoder::AdaptiveDecode_CRCAlternate(
    std::vector<double> received_signal,
    std::vector<std::chrono::milliseconds>& timeDurations) {
  /**
   * @brief
   *
   * @param received_signal
   * @param timeDurations
   * @return DLDInfo
   */

  std::vector<std::vector<MessageInformation>> output;
  std::vector<MessageInformation> output_0;
  std::vector<MessageInformation> output_1;
  DualListMap mp;
  double best_current_match = INT_MAX;
  double decoder_threshold_0 = INT_MAX;
  double decoder_threshold_1 = INT_MAX;

  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];

  // crc traceback degree comparison
  // the smaller crc should traceback once every
  int list_decoder_0_num_tracebacks = 1;
  int list_decoder_1_num_tracebacks = 1;
  if (code_0.v <= code_1.v) {
    list_decoder_1_num_tracebacks = crc_ratio_;
  } else {
    list_decoder_0_num_tracebacks = crc_ratio_;
  }

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

  assert(received_codec_1.size() == 74);
  assert(received_codec_2.size() == 74);

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      ConstructZTCCTrellis_WithList_ProductMetric(
          received_codec_1, code_0, trellis_ptrs_[0], timeDurations[0]);
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
      ConstructZTCCTrellis_WithList_ProductMetric(
          received_codec_2, code_1, trellis_ptrs_[1], timeDurations[0]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;
  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;

  // time taken to do additional traceback and insertion
  Stopwatch traceback_insertion_stopwatch;
  traceback_insertion_stopwatch.tic();
  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 >= max_path_to_search_) {
      decoder_0_stop = true;
    }
    if (!decoder_0_stop) {
      // std::cout << "decoder 0 traceback" << std::endl;
      std::vector<MessageInformation> mi_0_out(list_decoder_0_num_tracebacks);
      mi_0_out =
          TraceBack_Multiple(heap_list_0, code_0, trellis_ptrs_[0], trellis_0,
                             prev_paths_list_0, num_path_searched_0,
                             num_total_stages_0, list_decoder_0_num_tracebacks);
      for (MessageInformation mi_0 : mi_0_out) {
        mi_0.decoder_index = 0;
        if (!mi_0.list_size_exceeded) {
          // if list size is not exceeded
          if (mi_0.path_metric >= decoder_threshold_0) {
            decoder_0_stop = true;
          }
          mp.insert(mi_0);
          output_0.push_back(mi_0);
        }
      }
    }
    // list decoder 1 traceback
    if (num_path_searched_1 >= max_path_to_search_) {
      decoder_1_stop = true;
    }
    if (!decoder_1_stop) {
      // std::cout << "decoder 1 traceback" << std::endl;
      // std::cout << "before " << num_path_searched_1 << std::endl;
      std::vector<MessageInformation> mi_1_out(list_decoder_1_num_tracebacks);
      mi_1_out =
          TraceBack_Multiple(heap_list_1, code_1, trellis_ptrs_[1], trellis_1,
                             prev_paths_list_1, num_path_searched_1,
                             num_total_stages_1, list_decoder_1_num_tracebacks);
      // std::cout << "after " << num_path_searched_1 << std::endl;
      for (MessageInformation mi_1 : mi_1_out) {
        mi_1.decoder_index = 1;
        if (!mi_1.list_size_exceeded) {
          if (mi_1.path_metric >= decoder_threshold_1) {
            decoder_1_stop = true;
          }
          mp.insert(mi_1);
          output_1.push_back(mi_1);
        }
      }
    }

    if (decoder_0_stop && decoder_1_stop) {
      // std::cout << "Both decoders declared stop" << std::endl;
      break;
    }

    if (mp.queue_size() != 0) {
      // std::cout << "Agreed message found: " << std::endl;
      // std::cout << "CRC alternate algorithm ended with path searched: [" <<
      // num_path_searched_0 << ", " << num_path_searched_1 << std::endl;
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // record the received signal just in case the decoding is incorrect
      agreed_message.received_signal = received_signal;
      decoder_threshold_0 = best_current_match - node_1.path_metric;
      decoder_threshold_1 = best_current_match - node_0.path_metric;
      best_combined_found = true;
      // free pointers
      delete heap_list_0;
      delete heap_list_1;
      heap_list_0 = nullptr;
      heap_list_1 = nullptr;
      // std::cout << "found agreed message" << std::endl;
      // std::cout << "Debug: agreed message list size: ";
      // dualdecoderutils::print(agreed_message.list_ranks);
      // std::cout << std::endl;

      return agreed_message;
    }
  }
  traceback_insertion_stopwatch.toc();
  timeDurations[2] += traceback_insertion_stopwatch.getElapsed();
  traceback_insertion_stopwatch.reset();

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  // when the output is not found, we store the message to be all -1's
  // and the received_signal
  empty_message.message = std::vector<int>(64, -1);
  empty_message.received_signal = received_signal;

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  // std::cout << "found nothing!" << std::endl;
  return empty_message;
}

DLDInfo DualListDecoder::AdaptiveDecode_SimpleAlternate_rate_1_2(
    std::vector<double> received_signal,
    std::vector<std::chrono::milliseconds>& timeDurations) {
  /**
  @brief This function adaptively expand the list size until the smallest future
  match (SFM) metric is larger than the best current match (BCM) metric.

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

  // unleaver to unleave the bits from received_signal
  std::vector<double> received_codec_2;
  std::vector<double> received_codec_1;
  for (size_t i = 0; i < received_signal.size(); ++i) {
    if (i % 3 == 0) {
      received_codec_2.push_back(received_signal[i]);
      received_codec_2.push_back(received_signal[i + 1]);
    } else if (i % 3 == 1) {
      received_codec_1.push_back(received_signal[i]);
      received_codec_1.push_back(received_signal[i + 1]);
    } else {
      continue;
    }
  }
  assert(received_codec_2.size() == 148);
  assert(received_codec_1.size() == 148);
  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      ConstructZTCCTrellis_WithList_EuclideanMetric(
          received_codec_2, code_0, trellis_ptrs_[0], timeDurations[0]);
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
      ConstructZTCCTrellis_WithList_EuclideanMetric(
          received_codec_1, code_1, trellis_ptrs_[1], timeDurations[0]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;

  // time taken to do additional traceback and insertion
  Stopwatch sw_step2_3;
  sw_step2_3.tic();
  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 >= max_path_to_search_) {
      decoder_0_stop = true;
    }
    if (!decoder_0_stop) {
      // std::cout << "decoder 0 traceback" << std::endl;
      MessageInformation mi_0 = TraceBack_Single(
          heap_list_0, code_0, trellis_ptrs_[0], trellis_0, prev_paths_list_0,
          num_path_searched_0, num_total_stages_0);
      mi_0.decoder_index = 0;
      if (!mi_0.list_size_exceeded) {
        // if list size is not exceeded
        if (mi_0.path_metric >= decoder_threshold_0) {
          decoder_0_stop = true;
        }
        mp.insert(mi_0);
        output_0.push_back(mi_0);
      }
    }
    // list decoder 1 traceback
    if (num_path_searched_1 >= max_path_to_search_) {
      decoder_1_stop = true;
    }
    if (!decoder_1_stop) {
      // std::cout << "decoder 1 traceback" << std::endl;
      // std::cout << "before " << num_path_searched_1 << std::endl;
      MessageInformation mi_1 = TraceBack_Single(
          heap_list_1, code_1, trellis_ptrs_[1], trellis_1, prev_paths_list_1,
          num_path_searched_1, num_total_stages_1);
      // std::cout << "after " << num_path_searched_1 << std::endl;
      mi_1.decoder_index = 1;
      if (!mi_1.list_size_exceeded) {
        if (mi_1.path_metric >= decoder_threshold_1) {
          decoder_1_stop = true;
        }
        mp.insert(mi_1);
        output_1.push_back(mi_1);
      }
    }

    if (decoder_0_stop && decoder_1_stop) {
      // std::cout << "Both decoders declared stop" << std::endl;
      break;
    }

    if (mp.queue_size() != 0) {
      // std::cout << "Agreed message found: " << std::endl;
      // std::cout << "Simple alternate algorithm ended with path searched: ["
      // << num_path_searched_0 << ", " << num_path_searched_1 << std::endl;
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // record the received signal just in case the decoding is incorrect
      agreed_message.received_signal = received_signal;
      decoder_threshold_0 = best_current_match - node_1.path_metric;
      decoder_threshold_1 = best_current_match - node_0.path_metric;
      best_combined_found = true;
      // free pointers
      delete heap_list_0;
      delete heap_list_1;
      heap_list_0 = nullptr;
      heap_list_1 = nullptr;
      // std::cout << "found agreed message" << std::endl;
      // std::cout << "Debug: agreed message list size: ";
      // dualdecoderutils::print(agreed_message.list_ranks);
      // std::cout << std::endl;

      return agreed_message;
    }
    // std::cout << "hello" << std::endl;
  }
  sw_step2_3.toc();
  timeDurations[2] += sw_step2_3.getElapsed();
  sw_step2_3.reset();

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  // when the output is not found, we store the message to be all -1's
  // and the received_signal
  empty_message.message = std::vector<int>(64, -1);
  empty_message.received_signal = received_signal;
  // std::cout << "List size exceeded!" << std::endl;

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  // std::cout << "found nothing!" << std::endl;
  return empty_message;
}

DLDInfo DualListDecoder::LookAheadDecode_SimpleAlternate_rate_1_2(
    std::vector<double> received_signal,
    std::vector<std::chrono::milliseconds>& timeDurations,
    std::vector<double> metric_0, std::vector<double> metric_1) {
  /**
  @brief This function adaptively expand the list size until the smallest future
  match (SFM) metric is larger than the best current match (BCM) metric.

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

  // Hard decoding every symbol we receive (including puncturing)
  // get squares of the differences: a vector of length received_signal
  // compute min_add_zero = sum of all last of every three squares of
  // differences compute min_add_one = sum of all first of every three squares
  // of differences Once there is a match
  //      Compute full distance, d_hat by
  //                    1. reencoding the matched message usign rate 1/3 encoder
  //                    2. then subtract from the received vector
  //                    3. squares the differences and sum together all squares.

  //      Need to determine when to stop searching list 0
  //      stop searching list 0 when the next item in the list with list metric,
  //      d_0
  //            1. d_0 plus min_add_zero is greater to d_hat, then we stop
  //            searching on list 0.
  //               For anything not greater, compute crc... as if you havn't
  //               found a match
  //      Need to determine when to stop searching list 1
  //      stop searching list 1 when the next item in the list with list metric,
  //      d_1
  //            1. d_1 plus min_add_one is greater to d_hat, then we stop
  //            searching on list 1.
  //               For anything not greater, compute crc... as if you havn't
  //               found a match

  //      if a better match is found, reset d_hat. and continue the search on
  //      both lists.
  //   Result: unless we reach the end of the lists, we cannot fail if SSV
  //   succeeds.
  
  Stopwatch debug_sw;
  debug_sw.tic();
  // Hard Decode
  std::vector<double> hard_decoding_result = HardDecode(received_signal);

  // Compute the minimum possible squared difference
  std::vector<double> squared_differences =
      ComputeSquaredDifferences(hard_decoding_result, received_signal);

  // Once the latest metric on list decoder 0 plus min_add_zero is greater than
  // the best match, stop searching on list 0
  double min_add_zero = SumGroupIndexElements(squared_differences, 3, 2);
  double min_add_one = SumGroupIndexElements(squared_differences, 3, 0);

  // std::cout << "min_add_zero: " << min_add_zero << std::endl;
  // std::cout << "min_add_one: " << min_add_one << std::endl;

  std::vector<std::vector<MessageInformation>> output;
  std::vector<MessageInformation> output_0;
  std::vector<MessageInformation> output_1;
  DualListMap mp;
  double best_current_match = INT_MAX;
  double decoder_threshold_0 = INT_MAX;
  double decoder_threshold_1 = INT_MAX;

  // unleaver to unleave the bits from received_signal
  std::vector<double> received_codec_2;
  std::vector<double> received_codec_1;

  std::vector<double> helper_trellis_2;
  std::vector<double> helper_trellis_1;
  for (size_t i = 0; i < received_signal.size(); ++i) {
    if (i % 3 == 0) {
      received_codec_2.push_back(received_signal[i]);
      received_codec_2.push_back(received_signal[i + 1]);
      helper_trellis_2.push_back(0.0);
      helper_trellis_2.push_back(received_signal[i + 2]);
    } else if (i % 3 == 1) {
      received_codec_1.push_back(received_signal[i]);
      received_codec_1.push_back(received_signal[i + 1]);
      helper_trellis_1.push_back(received_signal[i - 1]);
      helper_trellis_1.push_back(0.0);
    } else {
      continue;
    }
  }
  assert(received_codec_2.size() == 148);
  assert(received_codec_1.size() == 148);
  assert(helper_trellis_1.size() == 148);
  assert(helper_trellis_2.size() == 148);

  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // Look Ahead variables
  double d_hat = INT_MAX;
  // skip reconstruction if mp.size hasn't changed
  int mp_size = 0;

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      ConstructZTCCTrellis_WithList_EuclideanMetric(
          received_codec_2, code_0, trellis_ptrs_[0],
          timeDurations[0]);  // using 1st and 2nd
  int num_total_stages_0 = trellis_0[0].size();
  std::vector<std::vector<int>> prev_paths_list_0;
  MinHeap* heap_list_0 = new MinHeap;
  DetourNode node_0;
  node_0.start_state = 0;
  node_0.path_metric = trellis_0[0][num_total_stages_0 - 1].pathMetric;
  heap_list_0->insert(node_0);
  int num_path_searched_0 = 0;
  bool decoder_0_stop = false;
  bool decoder_0_LSE = false;

  // // decoder the leftover bits and improve min_add_zero
  // std::vector<std::vector<Cell>> trellis_helper_0 =
  //     ConstructZTCCTrellis_WithList_EuclideanMetric(helper_trellis_2, code_0,
  //     trellis_ptrs_[0], timeDurations[0]);
  // int num_total_stages_helper_0 = trellis_helper_0[0].size();
  // std::vector<std::vector<int>> prev_paths_list_helper_0;
  // MinHeap* heap_list_helper_0 = new MinHeap;
  // DetourNode node_helper_0;
  // node_helper_0.start_state = 0;
  // node_helper_0.path_metric = trellis_helper_0[0][num_total_stages_helper_0 -
  // 1].pathMetric; heap_list_helper_0->insert(node_helper_0); int
  // num_path_searched_helper_0 = 0;

  // MessageInformation mi_helper_0 = TraceBack_Single(heap_list_helper_0,
  // code_0, trellis_ptrs_[0], trellis_0,
  //                   prev_paths_list_0, num_path_searched_0,
  //                   num_total_stages_0);

  // list decoder 1
  std::vector<std::vector<Cell>> trellis_1 =
      ConstructZTCCTrellis_WithList_EuclideanMetric(
          received_codec_1, code_1, trellis_ptrs_[1],
          timeDurations[0]);  // using 2nd and 3rd
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = new MinHeap;

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  bool decoder_1_stop = false;
  bool decoder_1_LSE = false;

  debug_sw.toc();
  timeDurations[3] += debug_sw.getElapsed();
  debug_sw.reset();

  // time taken to do additional traceback and insertion
  Stopwatch sw_step2_3;
  sw_step2_3.tic();
  while (!best_combined_found) {
    // list decoder 0 traceback
    if (num_path_searched_0 >= max_path_to_search_) {
      decoder_0_LSE = true;
    }
    if (!decoder_0_stop && !decoder_0_LSE) {
      // pass in a parameter into Traceback_single to stop further traceback
      MessageInformation mi_0 = TraceBack_Single(
          heap_list_0, code_0, trellis_ptrs_[0], trellis_0, prev_paths_list_0,
          num_path_searched_0, num_total_stages_0);
      mi_0.decoder_index = 0;
      if (!mi_0.list_size_exceeded) {
        // if list size is not exceeded
        if (mi_0.path_metric + min_add_zero > best_current_match ||
            mi_0.path_metric >= decoder_threshold_0) {
          // std::cout << "Decoder 0 declared stop after crc list size " <<
          // mi_0.list_rank << std::endl;
          decoder_0_stop = true;
        }
        // std::cout << "mi_0 path metric before true metric: " <<
        // mi_0.path_metric << std::endl; Re-encode process
        std::vector<int> reencode_message = mi_0.message;
        for (int i = 0; i < encoder_.v; ++i) {
          reencode_message.push_back(0);
        }
        std::vector<int> reencode_encoded_msg =
            encoder_trellis_ptr_->encode(reencode_message);
        // BPSK modulate
        std::vector<int> reencode_codeword =
            BPSK::modulate(reencode_encoded_msg);
        // compute d_hat
        d_hat = dualdecoderutils::euclideanDistance(received_signal,
                                                    reencode_codeword);
        mi_0.path_metric = d_hat;
        mp.insert(mi_0);
        best_current_match =
            (best_current_match > d_hat) ? d_hat : best_current_match;
        output_0.push_back(mi_0);
      }
    }
    // list decoder 1 traceback
    if (num_path_searched_1 >= max_path_to_search_) {
      decoder_1_LSE = true;
    }
    if (!decoder_1_stop && !decoder_1_LSE) {
      // pass in a parameter into Traceback_single to stop further traceback
      MessageInformation mi_1 = TraceBack_Single(
          heap_list_1, code_1, trellis_ptrs_[1], trellis_1, prev_paths_list_1,
          num_path_searched_1, num_total_stages_1);

      mi_1.decoder_index = 1;
      if (!mi_1.list_size_exceeded) {
        if (mi_1.path_metric + min_add_one > best_current_match ||
            mi_1.path_metric >= decoder_threshold_1) {
          decoder_1_stop = true;
        }
        // Re-encode process
        std::vector<int> reencode_message = mi_1.message;
        for (int i = 0; i < encoder_.v; ++i) {
          reencode_message.push_back(0);
        }
        std::vector<int> reencode_encoded_msg =
            encoder_trellis_ptr_->encode(reencode_message);
        // BPSK modulate
        std::vector<int> reencode_codeword =
            BPSK::modulate(reencode_encoded_msg);
        // compute d_hatx
        d_hat = dualdecoderutils::euclideanDistance(received_signal,
                                                    reencode_codeword);
        mi_1.path_metric = d_hat;
        mp.insert(mi_1);
        best_current_match =
            (best_current_match > d_hat) ? d_hat : best_current_match;
        output_1.push_back(mi_1);
      }
    }

    // if both decoder take too long to verify or find anything
    if (decoder_0_LSE || decoder_1_LSE) {
      // if list size exceeded during verification
      if (mp.queue_size() != 0) {
        DLDInfo attemp_message = mp.pop_queue();
        if (attemp_message.combined_metric == best_current_match) {
          // std::cout << "Champion metric: " << attemp_message.combined_metric
          // << std::endl; std::cout << "Champion list ranks: " <<
          // attemp_message.list_ranks[0] << ", " <<
          // attemp_message.list_ranks[1] << std::endl; free pointers
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return attemp_message;
        }
      }
      // if list size exceeded and found nothing
      // simply break
      break;
    }

    if (mp.queue_size() != 0) {
      // std::cout << "Agreed message found: " << std::endl;
      // std::cout << "Simple alternate algorithm ended with path searched: ["
      // << num_path_searched_0 << ", " << num_path_searched_1 << std::endl;

      if (mp_size != mp.queue_size()) {
        mp_size = mp.queue_size();
        // DLDInfo agreed_message = mp.get_top();
      }

      // both decoders have to declare stop of some kind
      if (decoder_0_LSE || decoder_1_LSE ||
          (decoder_0_stop && decoder_1_stop)) {
        DLDInfo agreed_message = mp.pop_queue();
        best_current_match = agreed_message.combined_metric;
        // record the received signal just in case the decoding is incorrect
        agreed_message.received_signal = received_signal;
        decoder_threshold_0 = best_current_match - node_1.path_metric;
        decoder_threshold_1 = best_current_match - node_0.path_metric;
        best_combined_found = true;
        // free pointers
        delete heap_list_0;
        delete heap_list_1;
        heap_list_0 = nullptr;
        heap_list_1 = nullptr;
        // return
        return agreed_message;
      } else {
        continue;
      }
    }
  }
  sw_step2_3.toc();
  timeDurations[2] += sw_step2_3.getElapsed();
  // std::cout << "Time Measured: " << timeDurations[2].count() << std::endl;
  sw_step2_3.reset();

  output.push_back(output_0);
  output.push_back(output_1);
  DLDInfo empty_message;
  empty_message.combined_metric = INT_MAX;
  empty_message.list_ranks = {max_path_to_search_, max_path_to_search_};
  // when the output is not found, we store the message to be all -1's
  // and the received_signal
  empty_message.message = std::vector<int>(64, -1);
  empty_message.received_signal = received_signal;

  // free pointers
  delete heap_list_0;
  delete heap_list_1;
  heap_list_0 = nullptr;
  heap_list_1 = nullptr;
  return empty_message;
}

MessageInformation DualListDecoder::TraceBack_Single(
    MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
    const std::vector<std::vector<Cell>>& trellis_states,
    std::vector<std::vector<int>>& prev_paths, int& num_path_searched,
    int num_total_stages) {
  /**
   * @brief Does a single traceback using pre-build trellis pointed to by
   * trellis_ptr
   *
   * @param heap:A pointer to a Minheap that stores all detour nodes, which
   * encapsulates path metric.
   *
   * @param code:A struct of CodeInformation that stores k, n, v, crc_degree,
   * and crc_length.
   *
   * @param trellis_ptr:trellis output mapping and nextStates mapping
   *
   * @param trellis_states:A pointer that points to a 2d vector of Cells that
   * store
   *
   * @param prev_paths
   * @param num_path_searched
   * @param num_total_stages
   * @return MessageInformation
   */

  MessageInformation mi;
  bool found_path = false;

  while (!found_path) {
    if (num_path_searched >= max_path_to_search_) {
      // declare list size exceeded
      mi.path_metric = INT_MAX;
      mi.list_rank = num_path_searched + 1;
      mi.list_size_exceeded = true;
      return mi;
    }
    DetourNode detour = heap->pop();
    // TODO: check here if diff + partial
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

    if (CRC_Check(messageWithoutTrailingZeros, code.crc_length, code.crc_dec) &&
        path.front() == path.back() && path.back() == 0) {
      mi.path = path;
      mi.path_metric = detour.path_metric;
      // convert path to message
      mi.message = message;
      mi.list_rank = num_path_searched + 1;
      found_path = true;
    }

    num_path_searched++;
  }
  return mi;
}

std::vector<MessageInformation> DualListDecoder::TraceBack_Multiple(
    MinHeap* heap, const CodeInformation& code, FeedForwardTrellis* trellis_ptr,
    const std::vector<std::vector<Cell>>& trellis_states,
    std::vector<std::vector<int>>& prev_paths, int& num_path_searched,
    int num_total_stages, int num_trace_back) {
  /**
   * @brief Does multiple tracebacks using pre-build trellis pointed to by
   * trellis_ptr
   *
   * @param heap:A pointer to a Minheap that stores all detour nodes, which
   * encapsulates path metric.
   *
   * @param code:A struct of CodeInformation that stores k, n, v, crc_degree,
   * and crc_length.
   *
   * @param trellis_ptr:trellis output mapping and nextStates mapping
   *
   * @param trellis_states:A pointer that points to a 2d vector of Cells that
   * store
   *
   * @param prev_paths
   * @param num_path_searched
   * @param num_total_stages
   * @param num_trace_back
   * @return MessageInformation
   */

  std::vector<MessageInformation> output;
  for (int i = 0; i < num_trace_back; i++) {
    MessageInformation mi;
    bool found_path = false;

    while (!found_path) {
      if (num_path_searched >= max_path_to_search_) {
        // declare list size exceeded
        mi.list_size_exceeded = true;
        output.push_back(mi);
        return output;
      }
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
        double prev_path_metric =
            trellis_states[cur_state][stage - 1].pathMetric;
        forward_partial_path_metric += cur_path_metric - prev_path_metric;
        path[stage - 1] = cur_state;
      }
      prev_paths.push_back(path);

      std::vector<int> full_message = convertPathtoMessage(path, trellis_ptr);
      std::vector<int> messageWithoutTrailingZeros =
          convertPathtoTrimmedMessage(path, code, trellis_ptr);
      std::vector<int> message =
          deconvolveCRC(messageWithoutTrailingZeros, code);

      if (CRC_Check(messageWithoutTrailingZeros, code.crc_length,
                    code.crc_dec) &&
          path.front() == path.back() && path.back() == 0) {
        mi.path = path;
        mi.path_metric = detour.path_metric;
        // convert path to message
        mi.message = message;
        mi.list_rank = num_path_searched + 1;
        found_path = true;
      }

      num_path_searched++;
    }
    output.push_back(mi);
  }
  return output;
}



std::vector<std::vector<Cell>> DualListDecoder::constructZTListTrellis_precompute(
      const std::vector<double>& received_signal, CodeInformation code,
      FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time) {
  std::vector<std::vector<Cell>> trellisInfo;
  int lowrate_pathLength = (received_signal.size() / code.n) + 1;
  int lowrate_numStates = std::pow(2, code.v);

  trellisInfo = std::vector<std::vector<Cell>>(
      lowrate_numStates, std::vector<Cell>(lowrate_pathLength));

  // initializes just the zeroth valid starting states
  trellisInfo[0][0].pathMetric = 0;
  trellisInfo[0][0].init = true;

  // precomputing euclidean distance between the received signal and +/-1
  // std::vector<std::vector<double>> precomputedMetrics;
  // precomputedMetrics = std::vector<std::vector<double>>(received_signal.size(), std::vector<double>(2));
  
  // for (int state = 0; state < lowrate_numStates; state++) {
  //   for (int forwardPathIndex = 0; forwardPathIndex < std::pow(code.k,2); forwardPathIndex++) {
      
  //   }
  // }
  // for (int stage = 0; stage < received_signal.size(); stage++) {
  //   precomputedMetrics[stage][0] = std::pow(received_signal[stage] - 1, 2); // distance between +1
  //   precomputedMetrics[stage][1] = std::pow(received_signal[stage] + 1, 2); // distance between -1
  // }

  Stopwatch sw;
  sw.tic();
  // building the trellis
  for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
    for (int currentState = 0; currentState < lowrate_numStates;
         currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init) continue;
      
      
      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0;
           forwardPathIndex < trellis_ptr->nextStates_[0].size();
           forwardPathIndex++) {
        
        // we will only constrain to travelling on 0 path in the ending v stages
        // if (stage >= (lowrate_pathLength - code.v - 1)) {
        //   if (forwardPathIndex == 1) continue;
        // }

        // since our transitions correspond to symbols, the forwardPathIndex has
        // no correlation beyond indexing the forward path
        
        int nextState =
            trellis_ptr->nextStates_[currentState][forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0) continue;
        
        
        // double totalPathMetric =
        //     precomputedMetrics[stage][forwardPathIndex] + trellisInfo[currentState][stage].pathMetric;
        // // if (stage == 0) {
        // //   std::cout << "Debug: totalPathMetric: " << totalPathMetric << std::endl;
        // // }

        double branchMetric = 0;
        // std::cout << "code.n: " << code.n << std::endl;
        std::vector<int> output_point = dualdecoderutils::get_point(
            trellis_ptr->output_[currentState][forwardPathIndex], code.n);

        for (int i = 0; i < code.n; i++) {
          //// euclidean metric 
          // branchMetric += std::pow(
              // received_signal[code.n * stage + i] - (double)output_point[i], 2);

          //// absolute metric
          // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
          // stage + i] - (double)output_point[i]);

          //// special metric for BPSK
          branchMetric += -(received_signal[code.n * stage + i] * (double)output_point[i]);

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
          trellisInfo[nextState][stage + 1].subPathMetric =
              trellisInfo[nextState][stage + 1].pathMetric;
          trellisInfo[nextState][stage + 1].subFatherState =
              trellisInfo[nextState][stage + 1].fatherState;
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
        } else {
          trellisInfo[nextState][stage + 1].subPathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].subFatherState = currentState;
        }
      }
    }
  }
  sw.toc();
  ssv_time += sw.getElapsed();
  sw.reset();
  return trellisInfo;
}

std::vector<std::vector<Cell>>
DualListDecoder::ConstructZTCCTrellis_WithList_EuclideanMetric(
    const std::vector<double>& received_signal, CodeInformation code,
    FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time) {
  /**
  @brief Construct a ZTCC Trellis measuring the time taken by trellis
  construction (for both lists, iteratively) using a regular euclidean metric.

  @param received_signal: a vector of double that records the received values
  (deviates from +/-1).

  @param code: a dictionary encompasses all information related to a
  convolutional code.
        - 'n': convolutional code output number of bits. Example: n=2 for a rate
  1/2 code.
        - 'k': convolutional code input number of bits. Example: k=1 for a rate
  1/2 code.
        - 'v': number of memory elements in convolutional code. Example: v=6 for
  a CC with generator polynomial: x^6+x^5+x^4+x+1.

  @param trellis_ptr: a pointer that points to a FeedForwardTrellis object with
  following member variables
        - 'nextStates_': a 2d vector of size [2^(code.v), 2] for a binary
  convolutional code with v memory elements. At position (i, j), it determines
  the next state in integer when the currect state is i and currect input is j.
        - 'output_': a 2d vector of size [2^(code.v), 2] for a binary
  convolutional code with v memery elements. At position (i, j), it determines
  the output as integer (0 or 1 for binary code) when the currect state is i and
  currect input is j. To convert it to BPSK modulation, use
  dualdecoderutils::get_point().

  @param ssv_time: a chrono milliseconds object passed in by reference to
  increment an outer variable that keep track of dual list decoder trellis
  construction time. It tracks the time of add-compare-select operation of all
  active nodes.

  @note Theoretical complexity: 2^(v+1)-2 + 1.5*(2^(v+1)-2)
  + 1.5*(k+m-v)*2^(v+1) Define 1 unit of complexity as the complexity required
  to perform one addition.
         - 2^(v+1)-2: addition complexity on the first v section of a ZTCC
  trellis. Use sum of geometric series.
         - 1.5*(2^(v+1)-2): add-compare-select complexity on the last v section
  of a ZTCC trellis. Use sum of geometric series.
         - 1.5*(k+m-v)*2^(v+1): add-compare-select complexity on the middle
  (k+m-v) section of a ZTCC trellis.

  @return trellisInfo: std::vector<std::vector<Cell>>
          Each 'Cell' object constains the following information:
          -   struct Cell {
                            bool init = false;
                            double pathMetric = 3000;
                            int fatherState = -1;
                            double subPathMetric = 3000;
                            int subFatherState = -1;
                          };

          - 'init': a boolean variable indicating whether the Cell should be
  considered during forward propagation and traceback. Example: For a ZTCC, the
  lower left triangle of the first v stages and the lower right triangle of the
  last v stages do not need to be initialized.

          - 'pathMetric': a double variable that stores the optimal path metric
  up to that Cell.

          - 'fatherState': a integer variable that stores the optimal father
  state in integer of that Cell.

          - 'subPathMetric': a double variabel that stores the suboptimal path
  metric up to that Cell.

          - 'subFatherState': a integer variable that stores the suboptimal path
  metric up to that Cell.

          Note: The 'subPathMetric' and 'subFatherState' enable list decoding,
  i.e. additional traceback.

  */
  std::vector<std::vector<Cell>> trellisInfo;
  int lowrate_pathLength = (received_signal.size() / code.n) + 1;
  int lowrate_numStates = std::pow(2, code.v);
  
  Stopwatch sw;
  sw.tic();
  // repeat to measure time
  for (int i = 0; i < 1; i++) {
    trellisInfo = std::vector<std::vector<Cell>>(
        lowrate_numStates, std::vector<Cell>(lowrate_pathLength));

    // initializes just the zeroth valid starting states
    trellisInfo[0][0].pathMetric = 0;
    trellisInfo[0][0].init = true;

    // building the trellis
    for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
      for (int currentState = 0; currentState < lowrate_numStates;
          currentState++) {
        // if the state / stage is invalid, we move on
        if (!trellisInfo[currentState][stage].init) continue;

        // otherwise, we compute the relevent information
        for (int forwardPathIndex = 0;
            forwardPathIndex < trellis_ptr->nextStates_[0].size();
            forwardPathIndex++) {
          // we will only constrain to travelling on 0 path in the ending v stages
          // if (stage >= (lowrate_pathLength - code.v - 1)) {
          //   if (forwardPathIndex == 1) continue;
          // }

          // since our transitions correspond to symbols, the forwardPathIndex has
          // no correlation beyond indexing the forward path

          int nextState =
              trellis_ptr->nextStates_[currentState][forwardPathIndex];

          // if the nextState is invalid, we move on
          if (nextState < 0) continue;

          double branchMetric = 0;
          
          std::vector<int> output_point = dualdecoderutils::get_point(
              trellis_ptr->output_[currentState][forwardPathIndex], code.n);

          for (int i = 0; i < code.n; i++) {
            branchMetric += std::pow(
                received_signal[code.n * stage + i] - (double)output_point[i], 2);
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
            trellisInfo[nextState][stage + 1].subPathMetric =
                trellisInfo[nextState][stage + 1].pathMetric;
            trellisInfo[nextState][stage + 1].subFatherState =
                trellisInfo[nextState][stage + 1].fatherState;
            trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
            trellisInfo[nextState][stage + 1].fatherState = currentState;
          } else {
            trellisInfo[nextState][stage + 1].subPathMetric = totalPathMetric;
            trellisInfo[nextState][stage + 1].subFatherState = currentState;
          }
        }
      }
    }
  }

  sw.toc();
  ssv_time += sw.getElapsed();
  sw.reset();
  return trellisInfo;
}

std::vector<std::vector<Cell>>
DualListDecoder::ConstructZTCCTrellis_WithList_ProductMetric(
    const std::vector<double>& received_signal, CodeInformation code,
    FeedForwardTrellis* trellis_ptr, std::chrono::milliseconds& ssv_time) {
  /**
  @brief Construct a ZTCC Trellis measuring the time taken by trellis
  construction (for both lists, iteratively) using a special metric shown by
  Bill Ryan. Instead of calculating euclidean distance between received point
  and +/- 1, we compute the product of these two values.

  @param received_signal: a vector of double that records the received values
  (deviates from +/-1).

  @param code: a dictionary encompasses all information related to a
  convolutional code.
        - 'n': convolutional code output number of bits. Example: n=2 for a rate
  1/2 code.
        - 'k': convolutional code input number of bits. Example: k=1 for a rate
  1/2 code.
        - 'v': number of memory elements in convolutional code. Example: v=6 for
  a CC with generator polynomial: x^6+x^5+x^4+x+1.

  @param trellis_ptr: a pointer that points to a FeedForwardTrellis object with
  following member variables
        - 'nextStates_': a 2d vector of size [2^(code.v), 2] for a binary
  convolutional code with v memory elements. At position (i, j), it determines
  the next state in integer when the currect state is i and currect input is j.
        - 'output_': a 2d vector of size [2^(code.v), 2] for a binary
  convolutional code with v memery elements. At position (i, j), it determines
  the output as integer (0 or 1 for binary code) when the currect state is i and
  currect input is j. To convert it to BPSK modulation, use
  dualdecoderutils::get_point().

  @param ssv_time: a chrono milliseconds object passed in by reference to
  increment an outer variable that keep track of dual list decoder trellis
  construction time. It tracks the time of add-compare-select operation of all
  active nodes.

  @note Theoretical complexity: 2^(v+1)-2 + 1.5*(2^(v+1)-2)
  + 1.5*(k+m-v)*2^(v+1) Define 1 unit of complexity as the complexity required
  to perform one addition.
         - 2^(v+1)-2: addition complexity on the first v section of a ZTCC
  trellis. Use sum of geometric series.
         - 1.5*(2^(v+1)-2): add-compare-select complexity on the last v section
  of a ZTCC trellis. Use sum of geometric series.
         - 1.5*(k+m-v)*2^(v+1): add-compare-select complexity on the middle
  (k+m-v) section of a ZTCC trellis.

  @return trellisInfo: std::vector<std::vector<Cell>>
          Each 'Cell' object constains the following information:
          -   struct Cell {
                            bool init = false;
                            double pathMetric = 3000;
                            int fatherState = -1;
                            double subPathMetric = 3000;
                            int subFatherState = -1;
                          };

          - 'init': a boolean variable indicating whether the Cell should be
  considered during forward propagation and traceback. Example: For a ZTCC, the
  lower left triangle of the first v stages and the lower right triangle of the
  last v stages do not need to be initialized.

          - 'pathMetric': a double variable that stores the optimal path metric
  up to that Cell.

          - 'fatherState': a integer variable that stores the optimal father
  state in integer of that Cell.

          - 'subPathMetric': a double variabel that stores the suboptimal path
  metric up to that Cell.

          - 'subFatherState': a integer variable that stores the suboptimal path
  metric up to that Cell.

          Note: The 'subPathMetric' and 'subFatherState' enable list decoding,
  i.e. additional traceback.

  */
  std::vector<std::vector<Cell>> trellisInfo;
  int lowrate_pathLength = (received_signal.size() / code.n) + 1;
  int lowrate_numStates = std::pow(2, code.v);

  trellisInfo = std::vector<std::vector<Cell>>(
      lowrate_numStates, std::vector<Cell>(lowrate_pathLength));

  // initializes just the zeroth valid starting states
  trellisInfo[0][0].pathMetric = 0;
  trellisInfo[0][0].init = true;

  // precomputing euclidean distance between the received signal and +/-1
  std::vector<std::vector<double>> precomputedMetrics;
  precomputedMetrics = std::vector<std::vector<double>>(received_signal.size(),
                                                        std::vector<double>(2));

  for (int state = 0; state < lowrate_numStates; state++) {
    for (int forwardPathIndex = 0; forwardPathIndex < std::pow(code.k, 2);
         forwardPathIndex++) {
    }
  }
  for (int stage = 0; stage < received_signal.size(); stage++) {
    precomputedMetrics[stage][0] =
        std::pow(received_signal[stage] - 1, 2);  // distance between +1
    precomputedMetrics[stage][1] =
        std::pow(received_signal[stage] + 1, 2);  // distance between -1
  }

  Stopwatch sw;
  sw.tic();
  // building the trellis
  for (int stage = 0; stage < lowrate_pathLength - 1; stage++) {
    for (int currentState = 0; currentState < lowrate_numStates;
         currentState++) {
      // if the state / stage is invalid, we move on
      if (!trellisInfo[currentState][stage].init) continue;
      
      
      // otherwise, we compute the relevent information
      for (int forwardPathIndex = 0;
           forwardPathIndex < trellis_ptr->nextStates_[0].size();
           forwardPathIndex++) {
        
        // we will only constrain to travelling on 0 path in the ending v stages
        // if (stage >= (lowrate_pathLength - code.v - 1)) {
        //   if (forwardPathIndex == 1) continue;
        // }

        // since our transitions correspond to symbols, the forwardPathIndex has
        // no correlation beyond indexing the forward path
        
        int nextState =
            trellis_ptr->nextStates_[currentState][forwardPathIndex];

        // if the nextState is invalid, we move on
        if (nextState < 0) continue;
        
        
        // double totalPathMetric =
        //     precomputedMetrics[stage][forwardPathIndex] + trellisInfo[currentState][stage].pathMetric;
        // // if (stage == 0) {
        // //   std::cout << "Debug: totalPathMetric: " << totalPathMetric << std::endl;
        // // }

        double branchMetric = 0;
        // std::cout << "code.n: " << code.n << std::endl;
        std::vector<int> output_point = dualdecoderutils::get_point(
            trellis_ptr->output_[currentState][forwardPathIndex], code.n);

        for (int i = 0; i < code.n; i++) {
          //// euclidean metric 
          // branchMetric += std::pow(
              // received_signal[code.n * stage + i] - (double)output_point[i], 2);

          //// absolute metric
          // branchMetric += std::abs(receivedMessage[lowrate_symbolLength *
          // stage + i] - (double)output_point[i]);

          //// special metric for BPSK
          branchMetric += -(received_signal[code.n * stage + i] * (double)output_point[i]);

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
          trellisInfo[nextState][stage + 1].subPathMetric =
              trellisInfo[nextState][stage + 1].pathMetric;
          trellisInfo[nextState][stage + 1].subFatherState =
              trellisInfo[nextState][stage + 1].fatherState;
          trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].fatherState = currentState;
        } else {
          trellisInfo[nextState][stage + 1].subPathMetric = totalPathMetric;
          trellisInfo[nextState][stage + 1].subFatherState = currentState;
        }
      }
    }
  }
  sw.toc();
  ssv_time += sw.getElapsed();
  sw.reset();
  return trellisInfo;
}

std::vector<double> DualListDecoder::HardDecode(
    const std::vector<double>& received_signal) {
  std::vector<double> resultVector;

  for (const double& value : received_signal) {
    if (value > 0) {
      resultVector.push_back(1.0);
    } else if (value < 0) {
      resultVector.push_back(-1.0);
    } else {
      // If value is 0, you can decide how to handle it (set to 0.0 or some
      // default value)
      resultVector.push_back(0.0);
    }
  }

  return resultVector;
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
  // std::vector<int> result(output.size(), 0);
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

bool DualListDecoder::CRC_Check(std::vector<int> input_data, int crc_bits_num,
                                int crc_dec) {
  std::vector<int> CRC;
  dualdecoderutils::dec_to_binary(crc_dec, CRC, crc_bits_num);

  for (int ii = 0; ii <= (int)input_data.size() - crc_bits_num; ii++) {
    if (input_data[ii] == 1) {
      // Note: transform doesn't include .end
      std::transform(input_data.begin() + ii,
                     input_data.begin() + (ii + crc_bits_num), CRC.begin(),
                     input_data.begin() + ii, dualdecoderutils::bin_sum);
    }
  }
  bool zeros = std::all_of(input_data.begin(), input_data.end(),
                           [](int i) { return i == 0; });
  return zeros;
}