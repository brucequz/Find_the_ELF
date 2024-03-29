#ifndef MINHEAP_H
#define MINHEAP_H

#include <deque>  // Include deque for the MinHeap implementation
#include <vector>

struct DetourNode {
  DetourNode() : original_state(-1) {};
  double path_metric;  // pathMetric
  int original_state;  // 
  int start_state;     // ending / traceback starting state
  int original_path = -1;
  double forward_path_metric;
  int detour_stage;
  bool operator>(const DetourNode& other) const {
    return this->path_metric > other.path_metric;
  }
  bool operator<(const DetourNode& other) const {
    return this->path_metric < other.path_metric;
  }
};

class MinHeap {
 public:
  MinHeap();
  void insert(DetourNode node);
  DetourNode pop();
  
  // TODO: write another funciton that checks the top of the heap instead of popping it
  // return metric of the top of the heap
  double checkTop();

  int size() { return heap_.size(); }

 private:
  std::deque<DetourNode> heap_;
  void heapify(int index);
  int leftChildIndex(int index) { return 2 * index + 1; }
  int rightChildIndex(int index) { return 2 * index + 2; }
  int parentIndex(int index) { return (index - 1) / 2; }
};

#endif
