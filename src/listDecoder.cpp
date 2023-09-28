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

  

  return output;
}