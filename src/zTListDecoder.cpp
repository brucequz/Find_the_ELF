#include "../include/feedForwardTrellis.h"
#include "../include/minHeap.h"
#include "../include/viterbiCodec.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <climits>

namespace {

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
	for (int i=0; i<n; i++){
		bin_output[i] = -2 * bin_output[i] + 1;
	}
	return bin_output;
}
} // namespace

MessageInformation ViterbiCodec::ztListDecoding(std::vector<double> receivedMessage){
	std::vector<std::vector<Cell>> trellisInfo;
	int lowrate_pathLength = (receivedMessage.size() / n_) + 1;

	trellisInfo = std::vector<std::vector<Cell>>(numStates_, std::vector<Cell>(lowrate_pathLength));

	// initializes all the valid starting states
	for(int i = 0; i < numStates_; i++){
		trellisInfo[i][0].pathMetric = 0;
		trellisInfo[i][0].init = true;
	}
	
	// building the trellis
	for(int stage = 0; stage < lowrate_pathLength - 1; stage++){
		for(int currentState = 0; currentState < numStates_; currentState++){
			// if the state / stage is invalid, we move on
			if(!trellisInfo[currentState][stage].init)
				continue;

			// otherwise, we compute the relevent information
			for(int forwardPathIndex = 0; forwardPathIndex < trellis_ptr_->nextStates_[0].size(); forwardPathIndex++){
				// since our transitions correspond to symbols, the forwardPathIndex has no correlation 
				// beyond indexing the forward path

				int nextState = trellis_ptr_->nextStates_[currentState][forwardPathIndex];
				
				// if the nextState is invalid, we move on
				if(nextState < 0)
					continue;
				
				double branchMetric = 0;
				std::vector<int> output_point = get_point(trellis_ptr_->output_[currentState][forwardPathIndex], n_);
				
				for(int i = 0; i < n_; i++){
					branchMetric += std::pow(receivedMessage[n_ * stage + i] - (double)output_point[i], 2);
					// branchMetric += std::abs(receivedMessage[lowrate_symbolLength * stage + i] - (double)output_point[i]);
				}
				double totalPathMetric = branchMetric + trellisInfo[currentState][stage].pathMetric;
				
				// dealing with cases of uninitialized states, when the transition becomes the optimal father state, and suboptimal father state, in order
				if(!trellisInfo[nextState][stage + 1].init){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].fatherState = currentState;
					trellisInfo[nextState][stage + 1].init = true;
				}
				else if(trellisInfo[nextState][stage + 1].pathMetric > totalPathMetric){
					trellisInfo[nextState][stage + 1].pathMetric = totalPathMetric;
					trellisInfo[nextState][stage + 1].fatherState = currentState;
				}
			}

		}
	}

  // perform the single traceback
	MessageInformation output;

  double minMetric = INT_MAX;
  int endingState = -1;

  if(trellisInfo[0][lowrate_pathLength - 1].pathMetric < minMetric){
    minMetric = trellisInfo[0][lowrate_pathLength - 1].pathMetric;
    endingState = 0;
  }

  std::vector<int> path(lowrate_pathLength);
  path[lowrate_pathLength - 1] = endingState;
  int currentState = endingState;

  // traceback
  for(int stage = lowrate_pathLength - 1; stage > 0; stage--){
      currentState = trellisInfo[currentState][stage].fatherState;
      path[stage - 1] = currentState;
  }

  // std::vector<int> codeword = convertPathtoTrimmedMessage(path);
  std::vector<int> message = convertPathtoTrimmedMessage(path);
  //std::cout<< "decoded message: ";
  //for (int i=0; i<message.size(); i++){
  //  std::cout<< message[i];
  //}
  //std::cout<< " end of message with length " << message.size() << std::endl;

  output.message = message;
  output.path = path;
  return output;
}

