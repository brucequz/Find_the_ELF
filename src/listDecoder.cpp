#include "../include/feedForwardTrellis.h"
#include "../include/minHeap.h"
#include "../include/viterbiCodec.h"

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
  // std::cout << "message size: " << message.size() << std::endl;
  // std::cout << "truncating: " << v_ << std::endl;
  message.resize(message.size() - v_);
  return message;
}