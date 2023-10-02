#ifndef MINMAXHEAP_H
#define MINMAXHEAP_H

#include <vector>

struct DetourNode {
  DetourNode() : original_state(-1) {};
  int path_metric;  // pathMetric
  int original_state;
  int start_state;
  bool operator>(const DetourNode& other) const {
    return this->path_metric > other.path_metric;
  }
  bool operator<(const DetourNode& other) const {
    return this->path_metric < other.path_metric;
  }
};

class MinMaxHeap {
 public:
  MinMaxHeap() {};
  MinMaxHeap(int size);
  void build(std::vector<DetourNode> input);
  void insert(DetourNode node);
  DetourNode pop();
  int size() { return heap_.size(); }

 private:
  std::vector<DetourNode> heap_;

  // helper functions
  bool hasChild(int index) {return 2 * index + 1 < this->size();}
  
  void pushDownMin(std::vector<int> input, int index);
  void push_down_max(std::vector<int> input, int index);
  int leftChildIndex(int index) { return 2 * index + 1; }
  int rightChildIndex(int index) { return 2 * index + 2; }
  int parentIndex(int index) { return (index - 1) / 2; }
  int grandParentIndex(int index);
  std::vector<int> grandChildrenIndices(int index);

};



#endif