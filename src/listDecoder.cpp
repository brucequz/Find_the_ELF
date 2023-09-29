#include "feedForwardTrellis.h"
#include "minHeap.h"
#include "viterbiCodec.h"

std::vector<MessageInformation> ViterbiCodec::listViterbiDecoding(
    const std::vector<int>& coded) {
  std::vector<std::vector<Cell>> trellis_states = constructTrellis(coded);
  std::vector<MessageInformation> output;

  int num_total_stages = trellis_states[0].size();
  MinHeap heap;  // Detour Tree

  // add all final stage nodes to the heap
  for (int i = 0; i < numStates_; ++i) {
    DetourNode node;
    node.start_state = i;
    node.path_metric = trellis_states[i][num_total_stages - 1].pathMetric;
    heap.insert(node);
  }

  DetourNode min_node = heap.pop();
  std::vector<int> path(num_total_stages);
  std::vector<int> message(k_ * (num_total_stages - 1), 0);

  int cur_state = min_node.start_state;

  for (int stage = num_total_stages - 1; stage >= 0; --stage) {
    int father_state = trellis_states[cur_state][stage].fatherState;
    path[stage] = cur_state;

    if (trellis_states[cur_state][stage].subFatherState != -1) {
      DetourNode node;
      // TODO: Complete list decoder
    }

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

  return output;
}