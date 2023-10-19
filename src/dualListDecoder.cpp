#include "../include/dualListDecoder.h"

#include "../include/feedForwardTrellis.h"
#include <vector>
#include <cmath>

#include <unordered_map>



DualListDecoder::DualListDecoder(std::vector<CodeInformation> code_info)
    : code_info_(code_info) {
  CodeInformation code_list_0 = code_info[0];
  CodeInformation code_list_1 = code_info[1];
  
  FeedForwardTrellis* trellis_ptr_0 = new FeedForwardTrellis(code_list_0);
  FeedForwardTrellis* trellis_ptr_1 = new FeedForwardTrellis(code_list_0);
  trellis_ptrs_.push_back(trellis_ptr_0);
  trellis_ptrs_.push_back(trellis_ptr_1);



}

DualListDecoder::~DualListDecoder() {
  
}


std::vector<MessageInformation> DualListDecoder::adaptiveDecode(std::vector<double> received_signal) {
  /*
  This function adaptively expand the list size until the smallest future match (SFM) metric is 
  larger than the bext current match (BCM) metric.
  
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
  
  // Set up variables
  bool best_combined_found = false;
  CodeInformation code_0 = code_info_[0];
  CodeInformation code_1 = code_info_[1];


  // list decoder 0
  std::vector<std::vector<Cell>> trellis_0 = constructZTCCTrellis(received_signal, code_0, trellis_ptrs_[0]);
  int num_total_stages_0 = trellis_0[0].size();
  std::vector<std::vector<int>> prev_paths_list_0;
  MinHeap* heap_list_0 = &min_heaps_[0];

  DetourNode node_0;
  node_0.start_state = 0;
  node_0.path_metric = trellis_0[0][num_total_stages_0 - 1].pathMetric;
  heap_list_0->insert(node_0);
  int num_path_searched_0 = 0;

  // list decoder 1
  std::vector<std::vector<Cell>> trellis_1 = constructZTCCTrellis(received_signal, code_1, trellis_ptrs_[1]);
  int num_total_stages_1 = trellis_1[0].size();
  std::vector<std::vector<int>> prev_paths_list_1;
  MinHeap* heap_list_1 = &min_heaps_[1];

  DetourNode node_1;
  node_1.start_state = 0;
  node_1.path_metric = trellis_1[0][num_total_stages_1 - 1].pathMetric;
  heap_list_1->insert(node_1);
  int num_path_searched_1 = 0;
  
  while (!best_combined_found) {
    // list decoder 0 traceback
    

  }
  



  
}

MessageInformation DualListDecoder::traceBack(MinHeap* heap, CodeInformation code, FeedForwardTrellis* trellis_ptr, std::vector<std::vector<Cell>> trellis_states, std::vector<std::vector<int>>& prev_paths, 
                             int& num_path_searched, int num_total_stages) {
  MessageInformation mi;
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
  std::vector<int> messageWithoutTrailingZeros = convertPathtoTrimmedMessage(path, code, trellis_ptr);
  std::vector<int> message = deconvolveCRC(messageWithoutTrailingZeros, code);
  
  if (checkCRC(messageWithoutTrailingZeros, code) && path.front() == path.back() && path.back() == 0) {
    mi.path = path;
    mi.path_metric = detour.path_metric;
    // convert path to message
    mi.message = full_message;
    mi.list_rank = num_path_searched;
  }

  num_path_searched++;
  return mi;

}

std::vector<std::vector<Cell>> DualListDecoder::constructZTCCTrellis(const std::vector<double>& received_signal, CodeInformation code, FeedForwardTrellis* trellis_ptr) {
  /*
  Construct output trellis for traceback later with euclidean distances
  as path metrics.
  */

  
  std::vector<std::vector<Cell>> trellis_states;
  int signal_length = received_signal.size();
  int number_of_stages = signal_length / code.n;

  trellis_states.resize(std::pow(2, code.v), std::vector<Cell>(number_of_stages + 1));
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
        std::vector<int> expected_signal = bpsk::modulate(expected_output);

        double branch_metric =
            dualdecoderutils::euclideanDistance(target_message, expected_signal);
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
    const std::vector<int> path, CodeInformation code, FeedForwardTrellis* trellis_ptr) {
  std::vector<int> message;
  message = convertPathtoMessage(path, trellis_ptr);
  message.resize(message.size() - code.v);
  return message;
}

std::vector<int> DualListDecoder::deconvolveCRC(const std::vector<int>& output, CodeInformation code) {
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

bool DualListDecoder::checkCRC(std::vector<int> demodulated, CodeInformation code) {
  // check crc by dividing the demodulated signal with crc poly
  std::vector<int> crc_bin = CRC::decToBin(code.crc_dec, code.crc_length);

  for (int ii = 0; ii <= (int)demodulated.size() - code.crc_length; ii++) {
    if (demodulated[ii] == 1) {
      // Note: transform doesn't include .end
      std::transform(demodulated.begin() + ii,
                     demodulated.begin() + (ii + code.crc_length), crc_bin.begin(),
                     demodulated.begin() + ii, CRC::binSum);
    }
  }
  bool all_zero = std::all_of(demodulated.begin(), demodulated.end(),
                              [](int i) { return i == 0; });
  return all_zero;
}