#include <cassert>

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


std::vector<MessageInformation> ViterbiCodec::ZTCCListDecoding_fullInformation_WithConstraint(
  const std::vector<double>& received_signal) {
  /**
   * @brief ZTCC list decoding with euclidean distance as metrics
   * 
   * @param received_signal a vector of double
   * @return std::vector<MessageInformation> 
   * @todo check trellis output and nextStates
   */
  
  std::vector<std::vector<Cell>> trellisInfo = ConstructZTCCTrellis_WithList_EuclideanMetric(received_signal);
  int pathLength = trellisInfo[0].size();
  std::vector<MessageInformation> output;
  MinHeap heap;
  std::vector<std::vector<int>> previousPaths;

  DetourNode detour;
  detour.start_state = 0;
  detour.path_metric = trellisInfo[0][pathLength-1].pathMetric;
  heap.insert(detour);

  int numPathSearched = 0;
  int crc_passing_path = 0;
  while(numPathSearched < list_size_){
		DetourNode detour = heap.pop();
		std::vector<int> path(pathLength);

		int newTracebackStage = pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.start_state;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.original_path != -1){
			forwardPartialPathMetric = detour.forward_path_metric;
			newTracebackStage = detour.detour_stage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.original_path];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].subPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].subFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].subPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].subFatherState != -1){
				DetourNode localDetour;
				localDetour.detour_stage = stage;
				localDetour.original_path = numPathSearched;
				localDetour.path_metric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forward_path_metric = forwardPartialPathMetric;
				localDetour.start_state = 0;
				heap.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].fatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> ztmessage = convertPathtoTrimmedMessage(path);


		// ztcc decoding requires only a crc check since the starting / ending states are fixed
		if(CRC_Check(ztmessage, crc_length_, crc_dec_)){
            std::vector<int> original_message = deconvolveCRC(ztmessage);
			assert(path[0] == 0);
			assert(path.back() == 0);
			MessageInformation mi;
			mi.message = original_message;
			mi.path = path;
			mi.path_metric = detour.path_metric;
			mi.crc_passing_rank = crc_passing_path+1;
			mi.list_rank = numPathSearched+1;
			if (crc_passing_path % 1 == 0){
				output.push_back(mi);
			}
			crc_passing_path++;
		}
		numPathSearched++;
	}
  return output;
}

std::vector<MessageInformation> ViterbiCodec::ZTCCListDecoding_fullInformation_NoConstraint(
  const std::vector<double>& received_signal) {
  /**
   * @brief ZTCC list decoding with euclidean distance as metrics
   * 
   * @param received_signal a vector of double
   * @return std::vector<MessageInformation> 
   * @todo check trellis output and nextStates
   */
  
  std::vector<std::vector<Cell>> trellisInfo = ConstructZTCCTrellis_WithList_EuclideanMetric(received_signal);
  int pathLength = trellisInfo[0].size();
  std::vector<MessageInformation> output;
  MinHeap heap;
  std::vector<std::vector<int>> previousPaths;

  DetourNode detour;
  detour.start_state = 0;
  detour.path_metric = trellisInfo[0][pathLength-1].pathMetric;
  heap.insert(detour);

  int numPathSearched = 0;
  while(numPathSearched < list_size_){
		DetourNode detour = heap.pop();
		std::vector<int> path(pathLength);

		int newTracebackStage = pathLength - 1;
		double forwardPartialPathMetric = 0;
		int currentState = detour.start_state;

		// if we are taking a detour from a previous path, we skip backwards to the point where we take the
		// detour from the previous path
		if(detour.original_path != -1){
			forwardPartialPathMetric = detour.forward_path_metric;
			newTracebackStage = detour.detour_stage;

			// while we only need to copy the path from the detour to the end, this simplifies things,
			// and we'll write over the earlier data in any case
			path = previousPaths[detour.original_path];
			currentState = path[newTracebackStage];

			double suboptimalPathMetric = trellisInfo[currentState][newTracebackStage].subPathMetric;

			currentState = trellisInfo[currentState][newTracebackStage].subFatherState;
			newTracebackStage--;
			
			double prevPathMetric = trellisInfo[currentState][newTracebackStage].pathMetric;

			forwardPartialPathMetric += suboptimalPathMetric - prevPathMetric;
			
		}
		path[newTracebackStage] = currentState;

		// actually tracing back
		for(int stage = newTracebackStage; stage > 0; stage--){
			double suboptimalPathMetric = trellisInfo[currentState][stage].subPathMetric;
			double currPathMetric = trellisInfo[currentState][stage].pathMetric;

			// if there is a detour we add to the detourTree
			if(trellisInfo[currentState][stage].subFatherState != -1){
				DetourNode localDetour;
				localDetour.detour_stage = stage;
				localDetour.original_path = numPathSearched;
				localDetour.path_metric = suboptimalPathMetric + forwardPartialPathMetric;
				localDetour.forward_path_metric = forwardPartialPathMetric;
				localDetour.start_state = 0;
				heap.insert(localDetour);
			}
			currentState = trellisInfo[currentState][stage].fatherState;
			double prevPathMetric = trellisInfo[currentState][stage - 1].pathMetric;
			forwardPartialPathMetric += currPathMetric - prevPathMetric;
			path[stage - 1] = currentState;
		}
		previousPaths.push_back(path);

		std::vector<int> ztmessage = convertPathtoTrimmedMessage(path);
		std::vector<int> message = convertPathtoMessage(path);
		// std::cout<< "decoded message: " << std::endl;
		// for (int i=0; i<message.size(); i++){
		// 	std::cout<< message[i];
		// }
		// std::cout<< "end of message" << std::endl;

		// ztcc decoding requires only a crc check since the starting / ending states are fixed
		MessageInformation mi;
		mi.path_metric = detour.path_metric;
		mi.message = ztmessage;
		mi.path = path;
		mi.list_rank = numPathSearched+1;
		output.push_back(mi);

		numPathSearched++;
	}
  return output;
}

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