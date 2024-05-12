#include "../include/dualListDecoder.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

#include "../include/dualListMap.h"
#include "../include/feedForwardTrellis.h"
#include "../include/stopWatch.h"
#include "../include/CONSTANTS.h"

namespace{
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

double sumVector(const std::vector<double>& vec) {
  return std::accumulate(vec.begin(), vec.end(), 0.0);
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
  assert(received_codec_2.size() == BLOCK_SIZE);
  assert(received_codec_1.size() == BLOCK_SIZE);
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
      std::cout << "found agreed message" << std::endl;
      std::cout << "Debug: agreed message list size: ";
      dualdecoderutils::print(agreed_message.list_ranks);
      std::cout << std::endl;

      return agreed_message;
    }
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
  std::cout << "List size exceeded!" << std::endl;
  std::cout << "Max list size searched: " << max_path_to_search_ << std::endl;

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
  DLDInfo unmatched_best;
  unmatched_best.combined_metric = INT_MAX;
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
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_0.message;
          unmatched_best.list_ranks = {mi_0.list_rank, max_path_to_search_};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_1.message;
          unmatched_best.list_ranks = {max_path_to_search_, mi_1.list_rank};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
        else {
          // return the code with the best metric so far
          // but obviously it is un-matched
          assert(unmatched_best.combined_metric != INT_MAX);

          // clean-up
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return unmatched_best;
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

DLDInfo DualListDecoder::LookAheadDecode_SimpleAlternate_StopOnceMatchFound(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations) {
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
  DLDInfo unmatched_best;
  unmatched_best.combined_metric = INT_MAX;
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
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_0.message;
          unmatched_best.list_ranks = {mi_0.list_rank, max_path_to_search_};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_1.message;
          unmatched_best.list_ranks = {max_path_to_search_, mi_1.list_rank};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
        else {
          // return the code with the best metric so far
          // but obviously it is un-matched
          assert(unmatched_best.combined_metric != INT_MAX);

          // clean-up
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return unmatched_best;
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
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      
      if (best_current_match > unmatched_best.combined_metric) {
        // std::cout << "Matched found, but it does not have the best metric" << std::endl;
        // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;
        // std::cout << "Agreed metric: " << best_current_match << std::endl;

        // free pointers
        delete heap_list_0;
        delete heap_list_1;
        heap_list_0 = nullptr;
        heap_list_1 = nullptr;
        return unmatched_best;
      }
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

DLDInfo DualListDecoder::LookAheadDecode_SimpleAlternate_StopOnceMatchFound_WithListSizeExceeded(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations) {
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
  assert(received_codec_2.size() == BLOCK_SIZE);
  assert(received_codec_1.size() == BLOCK_SIZE);
  assert(helper_trellis_1.size() == BLOCK_SIZE);
  assert(helper_trellis_2.size() == BLOCK_SIZE);

  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // Look Ahead variables
  double d_hat = INT_MAX;
  DLDInfo unmatched_best;
  unmatched_best.combined_metric = INT_MAX;
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

        // TODO: add the 

        // BPSK modulate
        std::vector<int> reencode_codeword =
            BPSK::modulate(reencode_encoded_msg);
        std::vector<double> reencode_squared_differences =
          ComputeSquaredDifferences(reencode_codeword, received_signal);

        double symbol1_metric = SumGroupIndexElements(reencode_squared_differences, 3, 0);
        double symbol2_metric = SumGroupIndexElements(reencode_squared_differences, 3, 1);
        double symbol3_metric = SumGroupIndexElements(reencode_squared_differences, 3, 2);
        auto symbol_metrics = std::make_tuple(symbol1_metric, symbol2_metric, symbol3_metric);
        mi_0.symbol_metrics = symbol_metrics;

        // compute d_hat
        d_hat = symbol1_metric + symbol2_metric + symbol3_metric;
        mi_0.path_metric = d_hat;
        mp.insert(mi_0);
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_0.message;
          unmatched_best.symbol_metrics = symbol_metrics;
          unmatched_best.list_ranks = {mi_0.list_rank, max_path_to_search_};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
        std::vector<double> reencode_squared_differences =
          ComputeSquaredDifferences(reencode_codeword, received_signal);

        double symbol1_metric = SumGroupIndexElements(reencode_squared_differences, 3, 0);
        double symbol2_metric = SumGroupIndexElements(reencode_squared_differences, 3, 1);
        double symbol3_metric = SumGroupIndexElements(reencode_squared_differences, 3, 2);
        auto symbol_metrics = std::make_tuple(symbol1_metric, symbol2_metric, symbol3_metric);
        mi_1.symbol_metrics = symbol_metrics;

        // compute d_hat
        d_hat = symbol1_metric + symbol2_metric + symbol3_metric;
        mi_1.path_metric = d_hat;
        mp.insert(mi_1);
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_1.message;
          unmatched_best.symbol_metrics = symbol_metrics;
          unmatched_best.list_ranks = {max_path_to_search_, mi_1.list_rank};
          unmatched_best.received_signal = received_signal;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
          // attemp_message.list_ranks[1] << std::endl;
          // free pointers
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return attemp_message;
        } 
        else {
          // return the code with the best metric so far
          // but obviously it is un-matched
          assert(unmatched_best.combined_metric != INT_MAX);
          // std::cout << "Return unmatched best" << std::endl;
          // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;

          // clean-up
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return unmatched_best;
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
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // std::cout << "Best metric: " << best_current_match << std::endl;
      
      if (best_current_match > unmatched_best.combined_metric) {
        // std::cout << "Matched found, but it does not have the best metric" << std::endl;
        // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;
        // std::cout << "Agreed metric: " << best_current_match << std::endl;

        // free pointers
        delete heap_list_0;
        delete heap_list_1;
        heap_list_0 = nullptr;
        heap_list_1 = nullptr;
        return unmatched_best;
      }
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
  output.push_back(output_0);
  output.push_back(output_1);
  if (unmatched_best.combined_metric != INT_MAX) {
    // std::cout << "No agreed message found, but unmatched best found" << std::endl;
    // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;
    // std::cout << "unmatched best list ranks: " << unmatched_best.list_ranks[0] << ", " << unmatched_best.list_ranks[1] << std::endl;
    delete heap_list_0;
    delete heap_list_1;
    heap_list_0 = nullptr;
    heap_list_1 = nullptr;
    return unmatched_best;
  }
  // std::cout << "No agreed message found, no unmatched best found" << std::endl;
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


DLDInfo DualListDecoder::LookAheadDecode_SimpleAlternate_StopOnceMatchFound_WithListSizeExceeded_HalfMetricOnSharedSymbols(std::vector<double> received_signal, std::vector<std::chrono::milliseconds>& timeDurations) {
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
  assert(received_codec_2.size() == BLOCK_SIZE);
  assert(received_codec_1.size() == BLOCK_SIZE);
  assert(helper_trellis_1.size() == BLOCK_SIZE);
  assert(helper_trellis_2.size() == BLOCK_SIZE);

  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];
  // step 1 start: SSV add-compare-select and initial traceback time for both
  // decoders auto ssv_start_time = std::chrono::steady_clock::now();

  // Look Ahead variables
  double d_hat = INT_MAX;
  DLDInfo unmatched_best;
  unmatched_best.combined_metric = INT_MAX;
  // skip reconstruction if mp.size hasn't changed
  int mp_size = 0;

  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 =
      ConstructZTCCTrellis_WithList_EuclideanMetric_HalfSharedMetric(
          received_codec_2, code_0, trellis_ptrs_[0],
          timeDurations[0], 1);  // using 1st and 2nd
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
      ConstructZTCCTrellis_WithList_EuclideanMetric_HalfSharedMetric(
          received_codec_1, code_1, trellis_ptrs_[1],
          timeDurations[0], 0);  // using 2nd and 3rd
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

        std::vector<double> reencode_squared_differences =
          ComputeSquaredDifferences(reencode_codeword, received_signal);

        double symbol1_metric = SumGroupIndexElements(reencode_squared_differences, 3, 0);
        double symbol2_metric = SumGroupIndexElements(reencode_squared_differences, 3, 1);
        double symbol3_metric = SumGroupIndexElements(reencode_squared_differences, 3, 2);
        auto symbol_metrics = std::make_tuple(symbol1_metric, symbol2_metric, symbol3_metric);
        mi_0.symbol_metrics = symbol_metrics;

        // compute d_hat
        d_hat = symbol1_metric + symbol2_metric + symbol3_metric;
        
        mi_0.path_metric = d_hat;
        mp.insert(mi_0);
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_0.message;
          unmatched_best.list_ranks = {mi_0.list_rank, max_path_to_search_};
          unmatched_best.received_signal = received_signal;
          unmatched_best.symbol_metrics = symbol_metrics;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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

      // std::cout << "mi_1 metric: " << mi_1.path_metric << std::endl;

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

        std::vector<double> reencode_squared_differences =
          ComputeSquaredDifferences(reencode_codeword, received_signal);

        double symbol1_metric = SumGroupIndexElements(reencode_squared_differences, 3, 0);
        double symbol2_metric = SumGroupIndexElements(reencode_squared_differences, 3, 1);
        double symbol3_metric = SumGroupIndexElements(reencode_squared_differences, 3, 2);
        auto symbol_metrics = std::make_tuple(symbol1_metric, symbol2_metric, symbol3_metric);
        mi_1.symbol_metrics = symbol_metrics;

        // compute d_hat
        d_hat = symbol1_metric + symbol2_metric + symbol3_metric;

        mi_1.path_metric = d_hat;
        mp.insert(mi_1);
        if (d_hat < best_current_match) {
          unmatched_best.combined_metric = d_hat;
          unmatched_best.message = mi_1.message;
          unmatched_best.list_ranks = {max_path_to_search_, mi_1.list_rank};
          unmatched_best.received_signal = received_signal;
          unmatched_best.symbol_metrics = symbol_metrics;
          best_current_match = d_hat;
        }
        // best_current_match = (best_current_match > d_hat) ? d_hat : best_current_match;
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
          // attemp_message.list_ranks[1] << std::endl;
          // free pointers
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return attemp_message;
        } else {
          // return the code with the best metric so far
          // but obviously it is un-matched
          assert(unmatched_best.combined_metric != INT_MAX);
          // std::cout << "Return unmatched best" << std::endl;
          // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;

          // clean-up
          delete heap_list_0;
          delete heap_list_1;
          heap_list_0 = nullptr;
          heap_list_1 = nullptr;
          return unmatched_best;
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
      DLDInfo agreed_message = mp.pop_queue();
      best_current_match = agreed_message.combined_metric;
      // std::cout << "Best metric: " << best_current_match << std::endl;
      
      if (best_current_match > unmatched_best.combined_metric) {
        // std::cout << "Matched found, but it does not have the best metric" << std::endl;
        // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;
        // std::cout << "Agreed metric: " << best_current_match << std::endl;

        // free pointers
        delete heap_list_0;
        delete heap_list_1;
        heap_list_0 = nullptr;
        heap_list_1 = nullptr;
        return unmatched_best;
      }
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
  output.push_back(output_0);
  output.push_back(output_1);
  if (unmatched_best.combined_metric != INT_MAX) {
    // std::cout << "No agreed message found, but unmatched best found" << std::endl;
    // std::cout << "Best metric: " << unmatched_best.combined_metric << std::endl;
    // std::cout << "unmatched best list ranks: " << unmatched_best.list_ranks[0] << ", " << unmatched_best.list_ranks[1] << std::endl;
    // std::cout << "unmatched best symbol metrics: " << std::get<0>(unmatched_best.symbol_metrics) << ", " << std::get<1>(unmatched_best.symbol_metrics) << ", " << std::get<2>(unmatched_best.symbol_metrics) << std::endl;
    delete heap_list_0;
    delete heap_list_1;
    heap_list_0 = nullptr;
    heap_list_1 = nullptr;
    return unmatched_best;
  }
  // std::cout << "No agreed message found, no unmatched best found" << std::endl;
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