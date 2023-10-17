#include "../include/feedForwardTrellis.h"
#include "../include/minHeap.h"
#include "../include/viterbiCodec.h"

std::vector<int> ViterbiCodec::convertPathtoMessage(
    const std::vector<int> path) {
  std::vector<int> message;
  for (int path_id = 0; path_id < path.size() - 1; path_id++) {
    for (int forward_path = 0;
         forward_path < trellis_ptr_->nextStates_[0].size(); forward_path++) {
      if (trellis_ptr_->nextStates_[path[path_id]][forward_path] ==
          path[path_id + 1]) {
        message.push_back(forward_path);
      }
    }
  }
  return message;
}

std::vector<int> ViterbiCodec::convertPathtoTrimmedMessage(
    const std::vector<int> path) {
  std::vector<int> message;
  message = convertPathtoMessage(path);
  message.resize(message.size() - v_);
  return message;
}

std::vector<MessageInformation> ViterbiCodec::unconstraintListViterbiDecoding(
    const std::vector<double>& received_signal) {
  std::vector<MessageInformation> output;
  std::vector<std::vector<Cell>> trellis_states =
      constructTrellis(received_signal);

  int num_total_stages = trellis_states[0].size();
  std::vector<std::vector<int>> prev_paths;
  MinHeap heap;  // Detour Tree

  // add all final stage nodes to the heap
  // we don't have to insert all the paths that start from the last stage
  // say there are too many stages, then maybe the top ${listSize} number of
  // paths should be recorded
  for (int i = 0; i < numStates_; ++i) {
    DetourNode node;
    node.start_state = i;
    node.path_metric = trellis_states[i][num_total_stages - 1].pathMetric;
    heap.insert(node);
  }

  int num_path_searched = 0;
  while (output.size() < list_size_) {
    DetourNode detour = heap.pop();
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
        localDetour.start_state = detour.start_state;
        localDetour.path_metric =
            cur_sub_path_metric + forward_partial_path_metric;
        localDetour.forward_path_metric = forward_partial_path_metric;
        localDetour.original_path = num_path_searched;
        localDetour.detour_stage = stage;
        heap.insert(localDetour);
      }

      cur_state = trellis_states[cur_state][stage].fatherState;
      double prev_path_metric = trellis_states[cur_state][stage - 1].pathMetric;
      forward_partial_path_metric += cur_path_metric - prev_path_metric;
      path[stage - 1] = cur_state;
    }
    prev_paths.push_back(path);
    std::vector<int> message = convertPathtoMessage(path);
    // message.resize(message.size() - v_);

    // if (checkCRC(message) && path.front() == path.back() && path.back() == 0) {
    //   MessageInformation mi;
    //   mi.path = path;
    //   mi.path_metric = detour.path_metric;
    //   // convert path to message
    //   mi.message = message;
    //   mi.list_rank = num_path_searched;
    //   output.push_back(mi);
    // }

    MessageInformation mi;
    mi.path = path;
    // convert path to message
    mi.message = message;
    mi.path_metric = detour.path_metric;
    mi.list_rank = num_path_searched;
    output.push_back(mi);

    num_path_searched++;
  }

  return output;
}

std::vector<MessageInformation> ViterbiCodec::listViterbiDecoding(
    const std::vector<double>& received_signal) {
  std::vector<MessageInformation> output;
  std::vector<std::vector<Cell>> trellis_states =
      constructTrellis(received_signal);

  int num_total_stages = trellis_states[0].size();
  std::vector<std::vector<int>> prev_paths;
  MinHeap heap;  // Detour Tree

  // add all final stage nodes to the heap
  // we don't have to insert all the paths that start from the last stage
  // say there are too many stages, then maybe the top ${listSize} number of
  // paths should be recorded
  for (int i = 0; i < numStates_; ++i) {
    DetourNode node;
    node.start_state = i;
    node.path_metric = trellis_states[i][num_total_stages - 1].pathMetric;
    heap.insert(node);
  }

  int num_path_searched = 0;
  while (output.size() < list_size_) {
    DetourNode detour = heap.pop();
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
        localDetour.start_state = detour.start_state;
        localDetour.path_metric =
            cur_sub_path_metric + forward_partial_path_metric;
        localDetour.forward_path_metric = forward_partial_path_metric;
        localDetour.original_path = num_path_searched;
        localDetour.detour_stage = stage;
        heap.insert(localDetour);
      }

      cur_state = trellis_states[cur_state][stage].fatherState;
      double prev_path_metric = trellis_states[cur_state][stage - 1].pathMetric;
      forward_partial_path_metric += cur_path_metric - prev_path_metric;
      path[stage - 1] = cur_state;
    }
    prev_paths.push_back(path);
    std::vector<int> message = convertPathtoMessage(path);

    // TODO: check what's in the first 310 elements
    message.resize(message.size() - v_);

    if (checkCRC(message) && path.front() == path.back() && path.back() == 0) {
      MessageInformation mi;
      mi.path = path;
      mi.path_metric = detour.path_metric;
      // convert path to message
      mi.message = message;
      mi.list_rank = num_path_searched;
      output.push_back(mi);
    }

    // MessageInformation mi;
    // mi.path = path;
    // // convert path to message
    // mi.message = message;
    // mi.path_metric = detour.path_metric;
    // mi.list_rank = num_path_searched;
    // output.push_back(mi);

    num_path_searched++;
  }

  return output;
}

std::vector<MessageInformation> ViterbiCodec::ZTCCListViterbiDecoding(
    const std::vector<double>& received_signal) {
  /*
  Decoding paths that only end in the zero-th state. 
  
  */
  
  std::vector<MessageInformation> output;
  std::vector<std::vector<Cell>> trellis_states =
      constructZTCCTrellis(received_signal);

  int num_total_stages = trellis_states[0].size();
  std::vector<std::vector<int>> prev_paths;
  MinHeap heap;  // Detour Tree
  
  DetourNode node;
  node.start_state = 0;
  node.path_metric = trellis_states[0][num_total_stages - 1].pathMetric;
  heap.insert(node);

  int num_path_searched = 0;

  while (num_path_searched < list_size_) {
    DetourNode detour = heap.pop();
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
        heap.insert(localDetour);
      }

      cur_state = trellis_states[cur_state][stage].fatherState;
      double prev_path_metric = trellis_states[cur_state][stage - 1].pathMetric;
      forward_partial_path_metric += cur_path_metric - prev_path_metric;
      path[stage - 1] = cur_state;
    }
    prev_paths.push_back(path);
    

    std::vector<int> full_message = convertPathtoMessage(path);
    std::vector<int> messageWithoutTrailingZeros = convertPathtoTrimmedMessage(path);
    std::vector<int> message = deconvolveCRC(messageWithoutTrailingZeros);
    

    // TODO: can I divide and check crc??
    if (checkCRC(messageWithoutTrailingZeros) && path.front() == path.back() && path.back() == 0) {
      MessageInformation mi;
      mi.path = path;
      mi.path_metric = detour.path_metric;
      // convert path to message
      mi.message = full_message;
      mi.list_rank = num_path_searched;
      output.push_back(mi);
    }

    num_path_searched++;
  }
  
  return output;

}

std::vector<MessageInformation> ViterbiCodec::unconstraintZTCCDecoding(
    const std::vector<double>& received_signal) {
  /*
  Decoding paths that only end in the zero-th state. 
  
  */
  
  std::vector<MessageInformation> output;
  std::vector<std::vector<Cell>> trellis_states =
      constructZTCCTrellis(received_signal);

  int num_total_stages = trellis_states[0].size();
  std::vector<std::vector<int>> prev_paths;
  MinHeap heap;  // Detour Tree
  
  DetourNode node;
  node.start_state = 0;
  node.path_metric = trellis_states[0][num_total_stages - 1].pathMetric;
  heap.insert(node);

  int num_path_searched = 0;

  while (num_path_searched < list_size_) {
    DetourNode detour = heap.pop();
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
        heap.insert(localDetour);
      }

      cur_state = trellis_states[cur_state][stage].fatherState;
      double prev_path_metric = trellis_states[cur_state][stage - 1].pathMetric;
      forward_partial_path_metric += cur_path_metric - prev_path_metric;
      path[stage - 1] = cur_state;
    }
    prev_paths.push_back(path);
    std::vector<int> message = convertPathtoMessage(path);
    message.resize(message.size() - v_);
    

    MessageInformation mi;
    mi.path = path;
    // convert path to message
    mi.message = message;
    mi.path_metric = detour.path_metric;
    mi.list_rank = num_path_searched;
    output.push_back(mi);

    num_path_searched++;
  }
  
  return output;

}